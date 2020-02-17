//
//  ProblemConfig.h
//  ErrorEstimation
//
//  Created by Philippe Devloo on 09/05/18.
//

#ifndef ProblemConfig_h
#define ProblemConfig_h

#include <set>
#include "TPZAnalyticSolution.h"

/// class to guide the error estimator
struct ProblemConfig
{
    
    /// geometric mesh on which the computational meshes are based
    TPZGeoMesh *gmesh = 0;
    /// polynomial order of the original mesh
    int porder = 1;
    /// increment in internal order of flux and pressure
    int orderp = 1;
    int orderq = 1;
    bool Iscontinuouspressure = true;
    int hdivmaismais =0;

    
    /// number of uniform refinements applied to the mesh
    int ndivisions = 1;
    
    int dimension = 0;
    
    STATE alpha=1;
    /// directory where the files will be stored
    std::string dir_name = ".";
    /// name identifying the problem
    std::string problemname;
    /// set of materialids in the mesh
    std::set<int> materialids;
    /// set of boundary condition material ids
    std::set<int> bcmaterialids;
    /// exact solution
    TLaplaceExample1 exact;


    ProblemConfig() {};

    ProblemConfig(const ProblemConfig &cp) : gmesh(cp.gmesh),
                                             porder(cp.porder),
                                        orderp(cp.orderp),
                                        orderq(cp.orderq),
                                              problemname(cp.problemname),
                                             materialids(cp.materialids),
                                             bcmaterialids(cp.bcmaterialids),
                                             exact(cp.exact),
                                             ndivisions(cp.ndivisions),
                                            Iscontinuouspressure(cp.Iscontinuouspressure),
                                            
                                             dimension(cp.dimension),hdivmaismais(cp.hdivmaismais),dir_name(cp.dir_name)
                                        
                                            
    {
    }

    ProblemConfig &operator=(const ProblemConfig &cp) {
        gmesh = cp.gmesh;
        porder = cp.porder;
        orderp = cp.orderp;
        orderq = cp.orderq;
        
        problemname = cp.problemname;
        materialids = cp.materialids;
        bcmaterialids = cp.bcmaterialids;
        exact = cp.exact;
        dimension = cp.dimension;
        Iscontinuouspressure = cp.Iscontinuouspressure;

        ndivisions = cp.ndivisions;
        hdivmaismais = cp.hdivmaismais;
        
        dir_name = cp.dir_name;
    
        return *this;
    }
};

#endif /* ProblemConfig_h */
