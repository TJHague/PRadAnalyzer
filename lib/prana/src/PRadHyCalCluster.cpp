//============================================================================//
// Abstract class for PRad HyCal Clustering                                   //
// Different reconstruction methods can be implemented accordingly            //
//                                                                            //
// Chao Peng, Weizhi Xiong                                                    //
// 09/28/2016                                                                 //
//============================================================================//

#include "PRadHyCalCluster.h"



PRadHyCalCluster::PRadHyCalCluster()
{
    // place holder
}

PRadHyCalCluster::~PRadHyCalCluster()
{
    // place holder
}

PRadHyCalCluster* PRadHyCalCluster::Clone()
const
{
    return new PRadHyCalCluster(*this);
}

void PRadHyCalCluster::FormCluster(std::vector<ModuleHit> &hs, std::vector<ModuleCluster> &cls)
const
{
    // to be implemented by methods
}

