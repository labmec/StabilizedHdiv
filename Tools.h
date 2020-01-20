//
//  Tools.hpp
//  StabilizedHdiv
//
//  Created by Denise De Siqueira on 20/01/20.
//

#ifndef Tools_hpp
#define Tools_hpp

#include <stdio.h>
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzbndcond.h"
#include "TPZInterfaceEl.h"
//
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
//
#include "pzanalysis.h"
//
#include "pzmultiphysicselement.h"
#include "pzmultiphysicscompel.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzbuildmultiphysicsmesh.h"

#include "pzpoisson3d.h"
#include "mixedpoisson.h"

#include "pzanalysis.h"

#include "TPZVTKGeoMesh.h"
//
#include "pzlog.h"
//
#include "pzhdivfull.h"
#include "pzaxestools.h"
#include "TPZCopySolve.h"
#include "pzstrmatrix.h"


#endif /* Tools_hpp */
//meio homogeneo
TPZGeoMesh *GMesh2(REAL Lx, REAL Ly,bool triang_elements);

//meio heterogeneo
TPZGeoMesh *GMesh3(bool triang_elements);

TPZCompMesh *CMeshFlux2(int pOrder, TPZGeoMesh *gmesh);
TPZCompMesh *CMeshPressure2(int pOrder,TPZGeoMesh *gmesh);
TPZCompMesh *CMeshMixed2(TPZVec<TPZCompMesh *> meshvec,TPZGeoMesh * gmesh);

void RefinamentoUnif(TPZGeoMesh* gmesh, int nDiv);
void ResolverSistema(TPZAnalysis &an, TPZCompMesh *fCmesh, int numthreads, bool direct);
void PosProcessMultph(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile);
void PosProcessFluxo(TPZAnalysis &an, std::string plotfile);

//solucao exata meio homogeneo
void SolExataMista(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux);
void SolExataPressao(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux);

//lado direito da equacao - meio homogeneo
void ForcingMista(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

//Para condicao de contorno - meio homogeneo
void NeumannAcima(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void NeumannAbaixo(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void DirichletEsquerda(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

//solucoes exata meio heterogeneo
void SolFluxoHeter(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux);
void SolPressaoHeter(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux);
void ForcingHeter(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

void NeumannAcimaXMaiorZero(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void NeumannAcimaXMenorZero(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void NeumannAbaixoXMaiorZero(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void NeumannAbaixoXMenorZero(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void DirichletXIgualMenosUm(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void DirichletXIgualUm(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void PermeabilityTensor(const TPZVec<REAL> &pt, TPZVec<STATE> &kabs, TPZFMatrix<STATE> &tensorK);



//erros
void ErrorHDiv2(TPZCompMesh *hdivmesh, std::ostream &out);
void ErrorL22(TPZCompMesh *l2mesh, std::ostream &out);

void ComputeFluxError(TPZCompMesh *cmesh, std::ostream &out);
void ComputePressureError(TPZCompMesh *cmesh, std::ostream &out);

TPZFMatrix<STATE> * ComputeInverse(TPZCompMesh * mphysics);

void NEquationsCondensed(TPZCompMesh *cmesh, long &neqglob,long &neqcond, bool ismisto);
TPZGeoMesh *CreateGeoMesh(int nel, TPZVec<int> &bcids);
