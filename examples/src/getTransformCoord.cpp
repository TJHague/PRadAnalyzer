//============================================================================//
// An example to search the best shift and tilting angles                     //
//                                                                            //
// Chao Peng                                                                  //
// 12/06/2017                                                                 //
//============================================================================//

#include "ConfigParser.h"
#include "generalstruct.h"
#include "canalib.h"
#include "PRadCoordSystem.h"
#include "PRadHyCalDetector.h"
#include <unordered_map>
#include <cmath>
#include <iostream>

using namespace std;

typedef Point3D<double> mypt;

#define HYCAL_Z (5732.27 - 88.9)
#define GEM1_Z (5314.55 - 88.9)
#define GEM2_Z (5274.91 - 88.9)
#define PROGRESS_COUNT 1000

static const mypt zaxis(0., 0., 1.);
static const mypt zero(0., 0., 0.);
static const mypt hycal_center(0., 0., HYCAL_Z);
// this value will be updated
static mypt gem_centers[2], gem_tilts[2];


ostream &operator <<(ostream &os, const mypt &p)
{
    os << p.x << ", " << p.y << ", " << p.z;
    return os;
}


struct PosDiff
{
    mypt diff, res, hycal, gem, gem1, gem2;
    unsigned int layout;
};

struct Ranger
{
    vector<double> min, max, step;
    vector<int> counts;
    int64_t size;
    size_t dim;

    Ranger(vector<double> min, vector<double> max, vector<int> counts)
    : min(min), max(max), counts(counts)
    {
        size = 1;
        dim = counts.size();
        for (size_t i = 0; i < dim; ++i)
        {
            size *= counts[i] + 1;
            if (counts[i] == 0) {
                step.push_back(0.);
            } else {
                step.push_back((max[i] - min[i])/(double)counts[i]);
            }
        }
        size += 1;
    }

    vector<double> At(int64_t count) const
    {
        vector<double> res;
        for(size_t i = 0; i < dim; ++i)
        {
            int curr = count%(counts[i] + 1);
            count /= counts[i] + 1;
            res.push_back(min[i] + step[i]*curr);
        }

        return res;
    }
};

mypt RevTransform(mypt trans, mypt rot, mypt coord);
unordered_map<string, PosDiff> ReadPositionDiffs(const string &path);
double TakeDifference(unordered_map<string, PosDiff> &diffs, const mypt &trans, const mypt &rot);


int main(int argc, char *argv[])
{
    // geometry
    auto obdiffs = ReadPositionDiffs("database/transform/delta_pos_table_2GeV.txt");
    unordered_map<string, PosDiff> gem[2];
    for (auto &it : obdiffs)
    {
        auto module = it.second.hycal;
        // edge at 44 mm, use 30 & 80 mm for cleaner events
        if (module.x < -80) {
            it.second.gem = it.second.gem1;
            gem[0][it.first] = it.second; // GEM1
        } else if (module.x > 80) {
            it.second.gem = it.second.gem2;
            gem[1][it.first] = it.second; // GEM2
        }
    }

    // search range
    // shift_x, shift_y, tilt_x, tilt_y, tilt_z
    vector<double> min = {-0.2, -0.2, 0, 0., -1.};
    vector<double> max = {0.2, 0.2, 0, 0., 1.};
    vector<int> counts= {20, 20, 0, 0, 100};

    Ranger range(min, max, counts);

    for (int k = 0; k < 2; ++k)
    {
        cout << "GEM " << k + 1 << endl;
        auto best_shift = gem_centers[k];
        auto best_tilt = gem_tilts[k];
        auto best_est = TakeDifference(gem[k], best_shift, best_tilt);
        auto best_diffs = gem[k];

        for (int64_t i = 0; i < range.size; ++i)
        {
            auto vals = range.At(i);
            auto shift = mypt(vals[0], vals[1], 0.) + gem_centers[k];
            auto tilt = mypt(vals[2], vals[3], vals[4]) + gem_tilts[k];
            auto est = TakeDifference(gem[k], shift, tilt);
            if(i%PROGRESS_COUNT == 0) {
                cout <<"------[ " << i << "/" << range.size << " ]---"
                     << "---[ estimator = "  << est << " ]---"
                     << "---[ shift and tilt = (" << shift << "), (" << tilt << ") ]------"
                     << "\r" << flush;
            }
            if (est < best_est) {
                best_shift = shift;
                best_tilt = tilt;
                best_est = est;
                best_diffs = gem[k];
            }
        }
        cout <<"------[ " << range.size << " points tested ]---"
             << "---[ best estimator = "  << best_est << " ]---"
             << "---[ shift and tilt = (" << best_shift << "), (" << best_tilt << ") ]------"
             << endl;
    }
     /*
    // test
    // point on hycal
    mypt hycal(20., 20., HYCAL_Z);
    // gem center position
    mypt gem1(1.5, 2.5, GEM1_Z);
    // gem tilt
    mypt gem1_tilt(100., 100., 100.);
    // normal of gem
    mypt norm = Transform(zero, gem1_tilt, zaxis);
    // intersection on gem plane of target->hycal projectile
    mypt p1 = hycal.intersect_plane(zero, gem1, norm);
    // gem coordinate
    mypt p2 = RevTransform(gem1, gem1_tilt, p1);
    // transform to general coordinate
    mypt p3 = Transform(gem1, gem1_tilt, p2);
    // project to hycal
    mypt p4 = p3.intersect_plane(zero, hycal, zaxis);

    cout << norm.x << ", " << norm.y << ", " << norm.z << endl;
    cout << p1.x << ", " << p1.y << ", " << p1.z << endl;
    cout << p2.x << ", " << p2.y << ", " << p2.z << endl;
    cout << p4.x << ", " << p4.y << ", " << p4.z << endl;
    */

    return 0;
}


// transform
mypt Transform(mypt trans, mypt rot, mypt coord)
{
    double x = coord.x;
    double y = coord.y;
    double z = coord.z;

    double xt, yt, zt;
    // firstly do the angle tilting
    // basic rotation matrix
    // Rx(a) = ( 1           0         0  )
    //         ( 0       cos(a)    sin(a) )
    //         ( 0      -sin(a)    cos(a) )
    xt = x, yt = y, zt = z;
    y = yt*cos(rot.x*0.001) + zt*sin(rot.x*0.001);
    z = -yt*sin(rot.x*0.001) + zt*cos(rot.x*0.001);

    // Ry(a) = ( cos(a)      0    -sin(a) )
    //         ( 0           1         0  )
    //         ( sin(a)      0     cos(a) )
    xt = x, yt = y, zt = z;
    x = xt*cos(rot.y*0.001) - zt*sin(rot.y*0.001);
    z = xt*sin(rot.y*0.001) + zt*cos(rot.y*0.001);

    // Rz(a) = ( cos(a)  sin(a)        0  )
    //         (-sin(a)  cos(a)        0  )
    //         ( 0           0         1  )
    xt = x, yt = y, zt = z;
    x = xt*cos(rot.z*0.001) + yt*sin(rot.z*0.001);
    y = -xt*sin(rot.z*0.001) + yt*cos(rot.z*0.001);

    // then correct the origin
    x += trans.x;
    y += trans.y;
    z += trans.z;

    return mypt(x, y, z);
}


// reversely transform
mypt RevTransform(mypt trans, mypt rot, mypt coord)
{
    double x = coord.x;
    double y = coord.y;
    double z = coord.z;

    trans = -1.*trans;
    rot = -1.*rot;

    // then correct the origin
    x += trans.x;
    y += trans.y;
    z += trans.z;

    double xt, yt, zt;
    // Rz(a) = ( cos(a)  sin(a)        0  )
    //         (-sin(a)  cos(a)        0  )
    //         ( 0           0         1  )
    xt = x, yt = y, zt = z;
    x = xt*cos(rot.z*0.001) + yt*sin(rot.z*0.001);
    y = -xt*sin(rot.z*0.001) + yt*cos(rot.z*0.001);

    // Ry(a) = ( cos(a)      0    -sin(a) )
    //         ( 0           1         0  )
    //         ( sin(a)      0     cos(a) )
    xt = x, yt = y, zt = z;
    x = xt*cos(rot.y*0.001) - zt*sin(rot.y*0.001);
    z = xt*sin(rot.y*0.001) + zt*cos(rot.y*0.001);

    // firstly do the angle tilting
    // basic rotation matrix
    // Rx(a) = ( 1           0         0  )
    //         ( 0       cos(a)    sin(a) )
    //         ( 0      -sin(a)    cos(a) )
    xt = x, yt = y, zt = z;
    y = yt*cos(rot.x*0.001) + zt*sin(rot.x*0.001);
    z = -yt*sin(rot.x*0.001) + zt*cos(rot.x*0.001);

    return mypt(x, y, z);
}


// read position differences
// convert differences to the module center coordinates on hycal, gem1, gem2
unordered_map<string, PosDiff> ReadPositionDiffs(const string &path)
{
    // read position difference data
    unordered_map<string, PosDiff> res;
    ConfigParser parser;
    parser.OpenFile(path);
    string name;
    PosDiff pdiff;
    while(parser.ParseLine())
    {
        parser >> name >> pdiff.diff.x >> pdiff.diff.y >> pdiff.res.x >> pdiff.res.y;
        res[name] = pdiff;
    }

    // transform difference back to detector coordinates
    // determine the difference between HyCal and GEMs
    PRadCoordSystem coord_sys("database/coordinates.dat");
    mypt trans1(0., 0., 0.), trans2(0., 0., 0.);
    mypt tilt1(0., 0., 0.), tilt2(0., 0., 0.);
    int count = 0;
    for (auto &coord : coord_sys.GetCoordsData())
    {
        auto det = coord.dets;
        count++;
        trans1 += det[PRadDetector::PRadGEM1].trans - det[PRadDetector::HyCal].trans;
        trans2 += det[PRadDetector::PRadGEM2].trans - det[PRadDetector::HyCal].trans;
        tilt1 += det[PRadDetector::PRadGEM1].rot;
        tilt2 += det[PRadDetector::PRadGEM2].rot;
    }
    trans1 /= (double)count;
    trans2 /= (double)count;
    tilt1 /= (double)count;
    tilt2 /= (double)count;
    // gem center
    gem_centers[0] = mypt(trans1.x, trans1.y, GEM1_Z);
    gem_tilts[0] = tilt1;
    cout << "GEM1 center: (" << gem_centers[0] << "), tilt: (" << gem_tilts[0] << ")" << endl;
    gem_centers[1] = mypt(trans2.x, trans2.y, GEM2_Z);
    gem_tilts[1] = tilt1;
    cout << "GEM2 center: (" << gem_centers[1] << "), tilt: (" << gem_tilts[1] << ")" << endl;

    // read hycal module geometry
    PRadHyCalDetector hycal;
    hycal.ReadModuleList("database/hycal_module.txt");
    // determine hycal and gem positions
    for (auto &p : res)
    {
        // hycal module center
        auto module = hycal.GetModule(p.first);
        p.second.hycal = mypt(module->GetX(), module->GetY(), module->GetZ());
        p.second.layout = module->GetLayoutFlag();

        // hycal module in the general coordinate system
        auto mcenter = Transform(hycal_center, zero, p.second.hycal);
        // project to hycal plane (lead glass module are not on plane)
        auto pr = mcenter.intersect_plane(zero, hycal_center, zaxis);
        // take difference (gem - hycal)
        auto gem_pos = pr + p.second.diff;
        // project back to gem plane, then convert to gem coordinates
        p.second.gem1 = RevTransform(gem_centers[0], gem_tilts[0], gem_pos.intersect_plane(zero, gem_centers[0], zaxis));
        p.second.gem2 = RevTransform(gem_centers[1], gem_tilts[1], gem_pos.intersect_plane(zero, gem_centers[1], zaxis));
    }
    return res;
}


double TakeDifference(unordered_map<string, PosDiff> &diffs,
                      const mypt &trans, const mypt &rot)
{
    double sigma = 0;
    for (auto it : diffs)
    {
        auto hycal = Transform(hycal_center, zero, it.second.hycal).intersect_plane(zero, hycal_center, zaxis);
        auto gem = Transform(trans, rot, it.second.gem).intersect_plane(zero, hycal_center, zaxis);
        auto diff = hycal - gem;
        sigma += cana::pow2(diff.x/it.second.res.x) + cana::pow2(diff.y/it.second.res.y);
    }

    return sqrt(sigma/(double)diffs.size());
}
