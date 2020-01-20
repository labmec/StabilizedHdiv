
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

#include "pzhdivfull.h"
#include "pzaxestools.h"
#include "TPZCopySolve.h"
#include "pzstrmatrix.h"

#include <iostream>
#include <math.h>

#include "Tools.h"
using namespace std;

//REAL const pi = 4.*atan(1.);
bool fTriang2 = false;
bool IsStab2 = true;
bool IsContinuou2 = true;
bool Useh22 = true;
//REAL Delta1 = 0.5;
//REAL Delta2 = 0.5;
bool IsFullHdiv2 = true;
bool IsHomogeneo2 = false;




//#ifdef LOG4CXX
//static LoggerPtr logger(Logger::getLogger("pz.hdiv"));
//#endif

#include "pztransfer.h"
int main(int argc, char *argv[])
{
#ifdef LOG4CXX
    std::string dirname = PZSOURCEDIR;
    std::string FileName = dirname;
    FileName = dirname + "/Projects/HDivEstabilizado/";
    FileName += "Hdivlog.cfg";
    InitializePZLOG(FileName);
#endif
    
    REAL Lx=0., Ly=0.;
    if(IsHomogeneo2==true){
        Lx = 1.;
        Ly = 0.5;
    }
    
    std::string outputfile("out");
    
    //ofstream saidaerro( "../erros-hdiv-estab.txt",ios::app);
    ofstream saidaerro( "erros-hdiv-estab.txt");
    
    saidaerro<<"\n";
    saidaerro <<"INFORMACOES DO TESTE: (0)->False e (1)->True \n";
    saidaerro <<"Malha Triangular: "<<fTriang2<<"\n";
    saidaerro <<"Problema eh homogeneo: "<< IsHomogeneo2<<"\n";
    saidaerro <<"Metodo estabilizado: " <<IsStab2<< ".    Usa parametro h2: "<< Useh22<<"\n";
    saidaerro <<"Espaco do fluxo eh full Hdiv: "<<IsFullHdiv2<<"\n";
    saidaerro <<"Espaco da pressao eh continuo: "<<IsContinuou2<<"\n\n";

    for(int p = 1; p<5; p++)
    {
        int pq = p;
        int pp;
        if(fTriang2==true){
            pp = pq-1;
        }else{
            pp = p;
        }
        
        int ndiv;
        saidaerro<<"\n CALCULO DO ERRO, COM ORDEM POLINOMIAL pq = " << pq << " e pp = "<< pp <<endl;
        for (ndiv = 1; ndiv < 6; ndiv++)
        {
            
            //std::cout << "p order " << p << " number of divisions " << ndiv << std::endl;
            
            saidaerro<<"\n<<<<<< Numero de divisoes uniforme ndiv = " << ndiv <<" >>>>>>>>>>> "<<endl;
            
            TPZGeoMesh *gmesh = NULL;
            if(IsHomogeneo2==true){
                gmesh = GMesh2(Lx, Ly,fTriang2);
            }else{
                gmesh = GMesh3(fTriang2);
            }
           // ofstream arg("gmesh1.txt");
            //gmesh->Print(arg);
            
            RefinamentoUnif(gmesh, ndiv);
            
            TPZCompMesh *cmesh1 = CMeshFlux2(pq,gmesh);
            TPZCompMesh *cmesh2 = CMeshPressure2(pp,gmesh);
//            
//            ofstream arg1("cmeshflux.txt");
//            cmesh1->Print(arg1);
//            
//            ofstream arg2("cmeshpressure.txt");
//            cmesh2->Print(arg2);

//            ofstream arg4("gmesh2.txt");
//            gmesh->Print(arg4);
            
            
            //malha multifisica
            TPZVec<TPZCompMesh *> meshvec(2);
            meshvec[0] = cmesh1;
            meshvec[1] = cmesh2;
            
            TPZCompMesh * mphysics = CMeshMixed2(meshvec,gmesh);
            
            mphysics->ExpandSolution();
            mphysics->CleanUpUnconnectedNodes();
            
            
//            ofstream arg4("gmeshMulti.txt");
//            mphysics->Print(arg4);
            
//            saidaerro << "Number of equations of flux " << cmesh1->NEquations() << std::endl;
//            saidaerro << "Number of equations of pressure " << cmesh2->NEquations() << std::endl;
//            saidaerro << "Number of equations TOTAL " << mphysics->NEquations() << "\n\n";
            
//            long neq_flux, neq_pres, neqcond_flux, neqcond_pres;
//            NEquationsCondensed(cmesh1, neq_flux, neqcond_flux);
//            NEquationsCondensed(cmesh2, neq_pres, neqcond_pres);
//            saidaerro << "Number of equations total flux: " <<neq_flux<< "\n";
//            saidaerro << "Number of equations condensadas flux: " <<neqcond_flux<< "\n";
//            saidaerro << "Number of equations total pressao: " <<neq_pres<< "\n";
//            saidaerro << "Number of equations condensadas pressao: " << neqcond_pres<< "\n\n";
            
            long neq_misto, neqcond_misto;
            NEquationsCondensed(mphysics, neq_misto, neqcond_misto,true);
            saidaerro << "Numero total de equacoes: " <<neq_misto<< "\n";
            saidaerro << "Numero de equacoes condensaveis: " <<neqcond_misto<< "\n";
            saidaerro << "Numero de equacoes final: " << neq_misto - neqcond_misto<< "\n\n";
/*
            int numthreads = 8;
            std::cout << "Number of threads " << numthreads << std::endl;

            
            bool opti = true;
            TPZAnalysis  * an = new TPZAnalysis(mphysics,opti);
            ResolverSistema(*an, mphysics,numthreads, true);
            
            
//            ofstream arg5("cmeshmultiphysics.txt");
//            mphysics->Print(arg5);

            //Calculo do erro
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
            TPZVec<REAL> erros;
    
            //saidaerro<<"\nErro da simulacao multifisica do fluxo (q)" <<endl;
            ErrorHDiv2(cmesh1, saidaerro);
            
            //saidaerro<<"\nErro da simulacao multifisica da pressao (p)" <<endl;
            ErrorL22(cmesh2, saidaerro);
            
            //Plot da solucao aproximada
            //string plotfile("Solution_mphysics.vtk");
            
//            char buf[256] ;
//            sprintf(buf,"ProblemaJuanGC_porder%d_h%d.vtk",p,ndiv);
//            PosProcessMultph(meshvec,  mphysics, *an, buf);
*/
        }
   
    }
    
	return EXIT_SUCCESS;
}

