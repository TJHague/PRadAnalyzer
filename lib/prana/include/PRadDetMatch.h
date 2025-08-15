//========================================================//
//match the HyCal Clusters and GEM clusters, both of them //
//need to be transformed to the Lab frame in advance      //
//========================================================//

#ifndef PRAD_DET_MATCH_H
#define PRAD_DET_MATCH_H

#include <string>
#include <vector>
#include "PRadEventStruct.h"
#include "ConfigObject.h"


class PRadDetMatch : public ConfigObject
{
public:
    PRadDetMatch(const std::string &path = "");
    virtual ~PRadDetMatch();

    void Configure(const std::string& path);

    std::vector<MatchHit> Match(std::vector<HyCalHit> &hycal,
                                const std::vector<GEMHit> &gem1,
                                const std::vector<GEMHit> &gem2) const;
    bool PreMatch(const HyCalHit &h, const GEMHit &g) const;
    void PostMatch(MatchHit &h) const;

    // Old matching algorithm causes an excess of fake matches
    //   in the GEM overlap region. E. Wrightson wrote an 
    //   algorithm that fixes this during analysis of sims.
    //   Merging into PRadAnalyzer - T.J. Hague Oct. 2024
    void PreMatchGEMs(const std::vector<GEMHit> &gem1,
                      const std::vector<GEMHit> &gem2) const;

private:
    float matchSigma;
    float overlapSigma;
    bool squareSel;
};

#endif
