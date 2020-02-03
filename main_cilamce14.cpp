
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

bool trapezoidalmesh = true;


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
    //ofstream saidaerro( "erros-hdiv-estab.txt",ios::app);
    
   // saidaerro<<"\n";
    
    ProblemConfig config;
    
    config.dimension = 2;
    TLaplaceExample1 example;
    //config.exact.fExact = example.EX;
    config.exact.fExact = example.ESinSin;
    config.Iscontinuouspressure = true;

    int maxp = 2;
    int maxhref = 5;

    TPZFMatrix<STATE> L2ErrorPressure(maxhref-1,maxp-1,0.);
    TPZFMatrix<STATE> EnergyErrorPressure(maxhref-1,maxp-1,0.);
    TPZFMatrix<STATE> L2ErrorFlux(maxhref-1,maxp-1,0.);
    TPZFMatrix<STATE> L2ErrorDivFlux(maxhref-1,maxp-1,0.);
    
    TPZFMatrix<STATE> L2ConvergPressure(maxhref-2,maxp-1,0.);
    TPZFMatrix<STATE> EnergyConvergPressure(maxhref-2,maxp-1,0.);
    TPZFMatrix<STATE> L2ConvergFlux(maxhref-2,maxp-1,0.);
    TPZFMatrix<STATE> L2ConvergDivFlux(maxhref-2,maxp-1,0.);
    
    config.vec_errors.Resize(4,0.);
    
    for(int p = 1; p< maxp; p++)
    {
        config.porder = p;
        config.orderp = p;
        config.orderq = p;


        for (int ndiv = 1; ndiv <maxhref; ndiv++)
        {
            config.ndivisions = ndiv;
            

            
            TPZGeoMesh *gmesh;
           
            if(trapezoidalmesh)
            {
                int nel= 2*ndiv;
                TPZManVector<int,4> bcids(4,-2);
                bcids[3] = -1;
                gmesh = CreateTrapezoidalMesh(nel, nel, 1.,1.,bcids);
                
                config.materialids.insert(1);
                config.bcmaterialids.insert(-1);//dirichlet
                config.bcmaterialids.insert(-2);//neumann
                config.gmesh = gmesh;
                gmesh->SetDimension(config.dimension);
                
                
                
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
                //bcids[1] = -1;
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
            
            L2ErrorPressure(ndiv-1,p-1) = config.vec_errors[0];
            EnergyErrorPressure(ndiv-1,p-1) = config.vec_errors[1];
            L2ErrorFlux(ndiv-1,p-1) = config.vec_errors[2];
            L2ErrorDivFlux(ndiv-1,p-1) = config.vec_errors[3];

        }

    }

    for(int j=0; j<maxp-1; j++){
        for(int i=0; i<maxhref-2; i++){
            L2ConvergPressure(i,j) = log(L2ErrorPressure(i+1,j)/L2ErrorPressure(i,j))/log(1./2.);
            EnergyConvergPressure(i,j) = log(EnergyErrorPressure(i+1,j)/EnergyErrorPressure(i,j))/log(1./2.);
            L2ConvergFlux(i,j) = log(L2ErrorFlux(i+1,j)/L2ErrorFlux(i,j))/log(1./2.);
            L2ConvergDivFlux(i,j) = log(L2ErrorDivFlux(i+1,j)/L2ErrorDivFlux(i,j))/log(1./2.);
        }
    }
    
    ofstream errtable;
    errtable.open("ErrosConvergencia.txt", ios::app);
    //std::ofstream errtable("ErrosConvergencia.txt",ios::app)
    errtable <<"\n\n";
    if(trapezoidalmesh){
    errtable <<"RESULTADOS PARA MALHA TRAPEZOIDAL"<<"\n\n";
    }else{
         errtable <<"RESULTADOS PARA MALHA QUADRILATERAL"<<"\n\n";
    }
    L2ErrorPressure.Print("Error L2 to pressure = ",errtable);
    L2ConvergPressure.Print("Convergence L2 to pressure = ",errtable);
    
    EnergyErrorPressure.Print("Error Energy to pressure  = ",errtable);
    EnergyConvergPressure.Print("Convergence Energy to pressure = ",errtable);
    
    L2ErrorFlux.Print("Error L2 to  flux = ",errtable);
    L2ConvergFlux.Print("Convergence L2 to  flux = ",errtable);
    
    L2ErrorDivFlux.Print("Error L2 to Divflux  = ",errtable);
    L2ConvergDivFlux.Print("Convergence L2 to Divflux  = ",errtable);
    errtable <<"--------------------------------------------------"<<"\n";
    errtable.close();
    
    return EXIT_SUCCESS;
}

