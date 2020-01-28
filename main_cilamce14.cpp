
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzbndcond.h"
#include "TPZInterfaceEl.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "TPZGeoLinear.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "pzgeoelside.h"

#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "TPZParSkylineStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "TPBSpStructMatrix.h"
#include "pzbstrmatrix.h"
#include "pzstepsolver.h"

#include "pzanalysis.h"

#include "pzmultiphysicselement.h"
#include "pzmultiphysicscompel.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzbuildmultiphysicsmesh.h"

#include "pzpoisson3d.h"
#include "mixedpoisson.h"

#include "pzanalysis.h"

#include "TPZVTKGeoMesh.h"

#include "pzlog.h"

//#include "pzhdivfull.h"
#include "pzaxestools.h"
#include "TPZCopySolve.h"
#include "pzstrmatrix.h"
#include "pzvtkmesh.h"

#include <iostream>
#include <math.h>

#include "Tools.h"
#include "ProblemConfig.h"
#include "TPZAnalyticSolution.h"

using namespace std;

//REAL const pi = 4.*atan(1.);
//bool fTriang2 = false;
//bool IsStab2 = true;
//bool IsContinuou2 = true;
//bool Useh22 = true;
//REAL Delta1 = 0.5;
//REAL Delta2 = 0.5;
//bool IsFullHdiv2 = true;
//bool IsHomogeneo2 = false;

bool trapezoidalmesh = false;


//#ifdef LOG4CXX
//static LoggerPtr logger(Logger::getLogger("pz.hdiv"));
//#endif

#include "pztransfer.h"
int main(int argc, char *argv[])
{
//#ifdef LOG4CXX
//    std::string dirname = PZSOURCEDIR;
//    std::string FileName = dirname;
//    FileName = dirname + "/Projects/HDivEstabilizado/";
//    FileName += "Hdivlog.cfg";
//    InitializePZLOG(FileName);
//#endif
    
    
    std::string outputfile("out");
    
    //ofstream saidaerro( "../erros-hdiv-estab.txt",ios::app);
    ofstream saidaerro( "erros-hdiv-estab.txt",ios::app);
    
    saidaerro<<"\n";
    
    ProblemConfig config;
    
    config.dimension = 2;
    TLaplaceExample1 example;
    config.exact.fExact = example.ESinSin;
    config.Iscontinuouspressure = true;
//    TPZManVector<int, 4> bcmaterialids(4,-1);
//    config.bcmaterialids.insert(-1);
//    config.bcmaterialids.insert(-2);
    

    
    for(int p = 1; p<2; p++)
    {
        config.porder = p;
        config.orderp = p;
        config.orderq = p;
        int pq = p;
        int pp = p;
        
        saidaerro<<"\n CALCULO DO ERRO, COM ORDEM POLINOMIAL pq = " << pq << " e pp = "<< pp <<endl;
<<<<<<< Updated upstream

=======
>>>>>>> Stashed changes
        for (int ndiv = 1; ndiv <2; ndiv++)
        {
            config.ndivisions = ndiv;
            
            //std::cout << "p order " << p << " number of divisions " << ndiv << std::endl;
            
            saidaerro<<"\n<<<<<< Numero de divisoes uniforme ndiv = " << ndiv <<" >>>>>>>>>>> "<<endl;
            
            TPZGeoMesh *gmesh;
           
            if(trapezoidalmesh)
            {
                gmesh = CreateTrapezoidalMesh(ndiv, ndiv, 1.,1.);
                
                ofstream arg("gmesh_trapez.txt");
                gmesh->Print(arg);
                
                std::ofstream file("GMeshTrapez.vtk");
                TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file);
            }
            else
            {
//                gmesh = CreateGeoMesh(1);
                TPZManVector<int,4> bcids(4,-2);
               // TPZManVector<int,4> bcids(4,-1);
                bcids[3] = -1;
                gmesh = CreateGeoMesh(1, bcids);
                config.materialids.insert(1);
                config.bcmaterialids.insert(-1);//dirichlet
                config.bcmaterialids.insert(-2);//neumann
                config.gmesh = gmesh;
                gmesh->SetDimension(config.dimension);
                ofstream arg("gmesh1.txt");
                gmesh->Print(arg);
                {
                    std::ofstream file("GMesh.vtk");
                    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file);
                }
                RefinamentoUnif(gmesh, ndiv);
                std::ofstream file("GMeshAposRef.vtk");
                TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file);
            }
            

            
            TPZMultiphysicsCompMesh *mphysics = CreateMultiphysicsMesh(config);
            
            mphysics->InitializeBlock();

            

            {
                std::ofstream out2("BuildComputationalMesh.vtk");
                TPZVTKGeoMesh::PrintCMeshVTK(mphysics, out2);
                ofstream arg4("BuildComputationalMesh.txt");
                mphysics->Print(arg4);
            }
            
            SolveStabilizedProblem(mphysics, config);
            

        }

    }

    return EXIT_SUCCESS;
}

