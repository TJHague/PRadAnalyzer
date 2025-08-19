//============================================================================//
// Match the HyCal Clusters and GEM clusters                                  //
// The clusters are required to be in the same frame (beam center frame)      //
// The frame transform work can be done in PRadCoordSystem                    //
//                                                                            //
// Details on projection:                                                     //
// By calling PRadCoordSystem::ProjectionDistance, the z from two clusters    //
// are compared first, if they are not at the same XOY-plane, the one with    //
// smaller z will be projected to the larger z. (from target point (0, 0, 0)) //
// If do not want the automatically projection, be sure you have projected    //
// the clusters to the same XOY-plane first.                                  //
//                                                                            //
//                                                                            //
// Xinzhan Bai, first version                                                 //
// Weizhi Xiong, adapted to PRad analysis software                            //
// Chao Peng, removed coordinates manipulation, only left matching part       //
// 10/22/2016                                                                 //
//============================================================================//

#include "PRadDetMatch.h"
#include "PRadCoordSystem.h"
#include <algorithm>
#include <set>

// constructor
PRadDetMatch::PRadDetMatch(const std::string &path)
{
    Configure(path);
}

// destructor
PRadDetMatch::~PRadDetMatch()
{
    // place holder
}

void PRadDetMatch::Configure(const std::string &path)
{
    // if no configuration file specified, load the default value quietly
    bool verbose = false;

    if(!path.empty()) {
        ConfigObject::Configure(path);
        verbose = true;
    }

    CONF_CONN(matchSigma, "Match_Factor", 5, verbose);
    CONF_CONN(overlapSigma, "GEM_Overlap_Factor", 10, verbose);
    CONF_CONN(squareSel, "Square_Selection", true, verbose);
}

template<typename T>
inline bool is_in(const std::set<T> &s, const T &i)
{
    return s.find(i) != s.end();
}

// define operator for set
// in this definition, if two gem hits are from the same gem detector and have
// the same reconstructed x, y, then they are the same hits
bool operator <(const GEMHit &lhs, const GEMHit &rhs)
{
    if(lhs.det_id != rhs.det_id)
        return lhs.det_id < rhs.det_id;

    if(lhs.x != rhs.x)
        return lhs.x < rhs.x;

    return lhs.y < rhs.y;
};

// Updated matching algorithm
// GEM planes are prematched to avoid double counting due to the overlap region
// The old prematch and postmatch functions are deprecated as keeping it in one function allows for substantial efficiency gains
std::vector<MatchHit> PRadDetMatch::Match(std::vector<HyCalHit>& hycal,
    const std::vector<GEMHit>& gem1,
    const std::vector<GEMHit>& gem2) const
{
    std::vector<MatchHit> results;
    // std::set<GEMHit> matched1, matched2;
    // std::vector<GEMHit> cand1, cand2;

    // sort in energy descending order
    std::sort(hycal.begin(), hycal.end(), [](const HyCalHit& h1, const HyCalHit& h2)
        {
            return h2.E < h1.E;
        });

    std::vector<GEMMatch> preMatchedGEMs = PreMatchGEMs(gem1, gem2);
    std::vector<bool> GEMMatchMask(preMatchedGEMs.size(), false);
    std::vector<std::tuple<GEMMatch, int, double>> GEMMatchDists; // GEMMatch, <index>, <distance>

    for (const auto& h : hycal)
    {
        GEMMatchDists.clear();

        // Old pre-matching logic
        float range = h.sig_pos * matchSigma;
        for (size_t i = 0; i < preMatchedGEMs.size(); i++)
        {
            if (GEMMatchMask[i]) continue;

            const GEMMatch& match = preMatchedGEMs[i];
            Point GEMs_At_HyCal = match.ProjectTo(h.z);
            double dist = GEMs_At_HyCal.dist(Point(h.x, h.y, h.z));

            if (dist < range) {
                GEMMatchDists.emplace_back(match, i, dist);
            }
        }

        // Post-matching logic
        // Find GEMMatchDists entry with shortest distance
        int best_index = -1;
        double best_dist = 1e9;
        for (int i = 0; i < GEMMatchDists.size(); i++) {
            if (std::get<2>(GEMMatchDists[i]) < best_dist) {
                best_dist = std::get<2>(GEMMatchDists[i]);
                best_index = i;
            }
        }

        // Create a matchhit object
        // If there is no matching GEM, create with an empty GEMMatch
        // If there is a matching GEM, create with the matched GEMs from all planes
        if(best_index < 0) {
            results.emplace_back(h, GEMMatch());
        } else {
            GEMMatchMask[std::get<1>(GEMMatchDists[best_index])] = true;
            results.emplace_back(h, GEMMatch(std::get<0>(GEMMatchDists[best_index])));
        }
    }

    return results;
}

// std::vector<MatchHit> PRadDetMatch::Match(std::vector<HyCalHit> &hycal,
//                                           const std::vector<GEMHit> &gem1,
//                                           const std::vector<GEMHit> &gem2)
// const
// {
//     std::vector<MatchHit> result;
//     std::set<GEMHit> matched1, matched2;
//     std::vector<GEMHit> cand1, cand2;

//     // sort in energy descendant order
//     std::sort(hycal.begin(), hycal.end(), [] (const HyCalHit &h1, const HyCalHit &h2)
//                                           {
//                                               return h2.E < h1.E;
//                                           });

//     for(size_t i = 0; i < hycal.size(); ++i)
//     {
//         const auto &hit = hycal.at(i);

//         // clean up cand1 and cand2 first
//         cand1.clear();
//         cand2.clear();

//         // pre match, only check if distance is within the range
//         // fill in hits as candidates
//         for(auto &ghit : gem1)
//         {
//             if(PreMatch(hit, ghit) && !is_in(matched1, ghit))
//                 cand1.push_back(ghit);
//         }
//         for(auto &ghit : gem2)
//         {
//             if(PreMatch(hit, ghit) && !is_in(matched2, ghit))
//                 cand2.push_back(ghit);
//         }

//         // no candidates
//         if(cand1.empty() && cand2.empty())
//             continue;

//         // create a new MatchHit
//         result.emplace_back(hit, std::move(cand1), std::move(cand2));
//         MatchHit &mhit = result.back();

//         // TODO remove it in PRadEventStruct.h
//         mhit.hycal_idx = i;

//         // find the matching status between HyCal and 2 GEMs
//         PostMatch(mhit);

//         // matched with gem1
//         if(TEST_BIT(mhit.mflag, kGEM1Match)) {
//             matched1.insert(mhit.gem1.front());
//         }
//         // matched with gem2
//         if(TEST_BIT(mhit.mflag, kGEM2Match)) {
//             matched2.insert(mhit.gem2.front());
//         }
//     }

//     return result;
// }

// project 1 HyCal cluster and 1 GEM cluster to HyCal Plane.
// if they are at the same z, there will be no projection, otherwise they are
// projected to the furthest z
// return true if they are within certain range (depends on HyCal resolution)
// return false if they are not
bool PRadDetMatch::PreMatch(const HyCalHit &hycal, const GEMHit &gem)
const
{
    // determine the matching range from hycal hit resolution and match factor
    float range = hycal.sig_pos*matchSigma;

    if(squareSel) {
        Point diff = PRadCoordSystem::ProjectionCoordDiff(hycal, gem);
        if(fabs(diff.x) > range || fabs(diff.y) > range)
            return false;
        else
            return true;
    } else {
        float dist = PRadCoordSystem::ProjectionDistance(hycal, gem);

        if(dist > range)
            return false;
        else
            return true;
    }
}

void PRadDetMatch::PostMatch(MatchHit &h)
const
{
    // sanity check, it should be rejected before calling the function
    if(h.gem1.empty() && h.gem2.empty()) {
        return;
    }

    // sort the candidates by delta_r between HyCal hit and GEM hits
    auto lamda_dr = [&h] (const GEMHit &h1, const GEMHit &h2)
                    {
                        return   PRadCoordSystem::ProjectionDistance(h.hycal, h1)
                               < PRadCoordSystem::ProjectionDistance(h.hycal, h2);
                    };

    std::sort(h.gem1.begin(), h.gem1.end(), lamda_dr);
    std::sort(h.gem2.begin(), h.gem2.end(), lamda_dr);

    // have candidates from gem1
    if(h.gem1.size()) {
        SET_BIT(h.mflag, kGEM1Match);
    }

    // have candidates from gem2
    if(h.gem2.size()) {
        SET_BIT(h.mflag, kGEM2Match);
    }

    // have candidates from both gem
    // need to check which one matches better
    if(h.gem1.size() && h.gem2.size()) {
        const GEMHit &hit1 = h.gem1.front(), &hit2 = h.gem2.front();
        float gem_dist = PRadCoordSystem::ProjectionDistance(hit1, hit2, PRadCoordSystem::target(), hit1.z);
        // not overlapping match
        if(gem_dist > overlapSigma * hit1.sig_pos) {
            float dist1 = PRadCoordSystem::ProjectionDistance(h, hit1);
            float dist2 = PRadCoordSystem::ProjectionDistance(h, hit2);
            // gem1 matches
            if(dist1 < dist2) {
                // remove the matching flag from gem2
                CLEAR_BIT(h.mflag, kGEM2Match);
            // gem2 matches
            } else {
                // remove the matching flag from gem1
                CLEAR_BIT(h.mflag, kGEM1Match);
            }
        }
    }

    // set the best matched gem hit and use gem position instead of hycal's
    // gem1 will always be used in overlapping match
    if(TEST_BIT(h.mflag, kGEM1Match)) {
        h.gem = h.gem1.front();
        h.SubstituteCoord(h.gem1.front());
    } else {
        h.gem = h.gem2.front();
        h.SubstituteCoord(h.gem2.front());
    }
}

std::vector<GEMMatch> PRadDetMatch::PreMatchGEMs(const std::vector<GEMHit>& gems1, const std::vector<GEMHit>& gems2) const {
    // If a gem is matched, set to true
    std::vector<bool> matchMask1(gems1.size(), false);
    std::vector<bool> matchMask2(gems2.size(), false);
    std::vector<GEMMatch> matchedGEMs;

    for(int i=0; i<gems1.size(); i++) {
        if(matchMask1[i]) continue;
        const GEMHit& g1 = gems1[i];
        // if(abs(GEM1_Layer1_Z - g1.z) > 1) continue; // Start with the first layer, get the second layer in loop 2
        // Check that it is in the first layer of the GEM

        // Storing best match details
        double best_dist = 1e9;
        int best_index = -1;

        for(int j=0; j<gems2.size(); j++) {
            if(matchMask2[j]) continue;
            const GEMHit& g2 = gems2[j];
            // if(abs(GEM1_Layer2_Z - g2.z) > 1) continue; // Check that it is in the second layer of the GEM

            // Implement the matching logic here
            double dist = PRadCoordSystem::ProjectionDistance(g1, g2, PRadCoordSystem::target(), g1.z);
            if(dist < matchSigma * g2.sig_pos && dist < best_dist) {
                best_dist = dist;
                best_index = j;
            }
        }

        if(best_index != -1) {
            matchMask1[i] = true;
            matchMask2[best_index] = true;

            // Create a matched GEM
            matchedGEMs.emplace_back(gems1[i], gems2[best_index]);
        }
        else {
            // Create a matched GEM with only one hit
            matchedGEMs.emplace_back(gems1[i]);
        }
    }

    return matchedGEMs;
}