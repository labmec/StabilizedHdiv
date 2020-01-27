
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

    
    for(int p = 1; p<2; p++)
    {
        config.porder = p;
        config.orderp = p;
        config.orderq = p;
        int pq = p;
        int pp = p;
        
        saidaerro<<"\n CALCULO DO ERRO, COM ORDEM POLINOMIAL pq = " << pq << " e pp = "<< pp <<endl;
        for (int ndiv = 2; ndiv < 5; ndiv++)
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
                gmesh = CreateGeoMesh(1);
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
            
            config.gmesh = gmesh;
            config.materialids.insert(1);
            config.bcmaterialids.insert(-1);
            gmesh->SetDimension(config.dimension);
            
//            TPZCompMesh *cmesh1 = CMeshFlux(config);//CMeshFlux(pq,gmesh);
//            ofstream arg1("cmesh_flux.txt");
//            cmesh1->Print(arg1);
//
//            TPZCompMesh *cmesh2 = CMeshPressure(config);//(pp,gmesh);

//            ofstream arg2("cmesh_pressure.txt");
//            cmesh2->Print(arg2);

            //malha multifisica
//            TPZVec<TPZCompMesh *> meshvec(2);
//            meshvec[0] = cmesh1;
//            meshvec[1] = cmesh2;
//            TPZCompMesh * mphysics = CMeshMixed(meshvec,config);//CMeshMixed(meshvec,gmesh);
            
            TPZMultiphysicsCompMesh *mphysics = CreateMultiphysicsMesh(config);
            
              mphysics->InitializeBlock();
            
            SolveStabilizedProblem(mphysics, config);
            

//            mphysics->ExpandSolution();
//            mphysics->CleanUpUnconnectedNodes();
            
//            ofstream arg4("cmesh_mixed.txt");
//            mphysics->Print(arg4);

//            saidaerro << "Number of equations of flux " << cmesh1->NEquations() << std::endl;
//            saidaerro << "Number of equations of pressure " << cmesh2->NEquations() << std::endl;
//            saidaerro << "Number of equations TOTAL " << cmesh1->NEquations()+cmesh2->NEquations()<<std::endl;
//            saidaerro << "Number of equations CONDENSADAS " << mphysics->NEquations() << "\n\n";

//            long neq_flux, neq_pres, neqcond_flux, neqcond_pres;
//           // NEquationsCondensed(cmesh1, neq_flux, neqcond_flux);
////            NEquationsCondensed(cmesh2, neq_pres, neqcond_pres);
//            saidaerro << "Number of equations total flux: " <<neq_flux<< "\n";
//            saidaerro << "Number of equations condensadas flux: " <<neqcond_flux<< "\n";
//            saidaerro << "Number of equations total pressao: " <<neq_pres<< "\n";
//            saidaerro << "Number of equations condensadas pressao: " << neqcond_pres<< "\n\n";
//
//            long neq_misto, neqcond_misto;
//            NEquationsCondensed(mphysics, neq_misto, neqcond_misto,true);
//            saidaerro << "Numero total de equacoes: " <<neq_misto<< "\n";
//            saidaerro << "Numero de equacoes condensaveis: " <<neqcond_misto<< "\n";
//            saidaerro << "Numero de equacoes final: " << neq_misto - neqcond_misto<< "\n\n";

//            int numthreads = 8;
//            std::cout << "Number of threads " << numthreads << std::endl;
//
//
//            bool opti = true;
//            TPZAnalysis  * an = new TPZAnalysis(mphysics,opti);
//            ResolverSistema(*an, mphysics,numthreads, true);


//            ofstream arg5("cmeshmultiphysics.txt");
//            mphysics->Print(arg5);

            //Calculo do erro
//            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
//            TPZVec<REAL> erros;
//
//            //saidaerro<<"\nErro da simulacao multifisica do fluxo (q)" <<endl;
//            ErrorHDiv2(cmesh1, saidaerro,config);
//
//            //saidaerro<<"\nErro da simulacao multifisica da pressao (p)" <<endl;
//            ErrorL22(cmesh2, saidaerro,config);

            //Plot da solucao aproximada
         //   string plotfile("Solution_mphysics.vtk");

//            char buf[256] ;
//            sprintf(buf,"ProblemaJuanGC_orderp%d_orderq_%d_h%d.vtk",config.orderp,config.orderq,config.ndivisions);
//            PosProcessMultph(meshvec,  mphysics, *an, buf);
//
        }

    }

    return EXIT_SUCCESS;
}

