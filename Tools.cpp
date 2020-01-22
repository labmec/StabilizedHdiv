//
//  Tools.cpp
//  StabilizedHdiv
//
//  Created by Denise De Siqueira on 20/01/20.
//

#include "Tools.h"
#include "pzgengrid.h"
#include "TPZAnalyticSolution.h"
#include "TPZMixedStabilizedHdiv.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.hdiv"));
#endif


int MatId =1;

int bcdirichlet = 0;
int bcneumann = 1;

int BC0=-1;
int BC1=-2;
int BC2=-3;
int BC3=-4;
int BC4=-5;
int BC5=-6;

REAL const pi = 4.*atan(1.);
bool fTriang = false;
bool IsStab = true;
bool IsContinuou = true;
bool Useh2 = false;
REAL Delta1 = 0.5;
REAL Delta2 = 0.5;
bool IsHomogeneo = true;

//bool IsFullHdiv = true;
//TRABALHO COM:
//static int bilinearounao [18] =   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//EM tpzquadrilateral.cpp de Topology

TPZGeoMesh *CreateGeoMesh(int nel) {
    
    TPZManVector<int> nx(2,nel);
    TPZManVector<REAL> x0(3,0.),x1(3,1.);
    x1[2] = 0.;
    TPZGenGrid gen(nx,x0,x1);
    gen.SetRefpatternElements(true);
    TPZGeoMesh* gmesh = new TPZGeoMesh;
    gen.Read(gmesh);
    gen.SetBC(gmesh, 4, BC0);
    gen.SetBC(gmesh, 5, BC1);
    gen.SetBC(gmesh, 6, BC2);
    gen.SetBC(gmesh, 7, BC3);
    
    //  UniformRefinement(1, gmesh);
    return gmesh;
}

TPZGeoMesh *CreateTrapezoidalMesh(int nelx, int nely, REAL Lx, REAL Ly){
    
    
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    
    TPZManVector<REAL,3> x0(3,0.),x1(3,0.);
    TPZManVector<int,3> nx(2);
    nx[0] = nelx;
    nx[1] = nely;
    x1[0] = Lx;
    x1[1] = Ly;
    TPZGenGrid gengrid(nx,x0,x1);
    
    gengrid.SetDistortion(0.25);
    //        gengrid.SetZigZagPattern();

    gengrid.Read(gmesh);
    x1[0] = Lx;
    x1[1] = 0.;
    gengrid.SetBC(gmesh, x0, x1, BC0);
    x0 = x1;
    x1[0] = Lx;
    x1[1] = Ly;
    gengrid.SetBC(gmesh, x0, x1, BC1);
    x0 = x1;
    x1[0] = 0.;
    x1[1] = Ly;
    gengrid.SetBC(gmesh, x0, x1, BC2);
    x0 = x1;
    x1[0] = 0.;
    x1[1] = 0.;
    gengrid.SetBC(gmesh, x0, x1, BC3);
    
    return gmesh;
}


TPZGeoMesh *GMesh2(REAL Lx, REAL Ly,bool triang_elements){
    
    int Qnodes = 4;
    
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    gmesh->SetMaxNodeId(Qnodes-1);
    gmesh->NodeVec().Resize(Qnodes);
    TPZVec<TPZGeoNode> Node(Qnodes);
    
    TPZVec <long> TopolQuad(4);
    TPZVec <long> TopolTriang(3);
    TPZVec <long> TopolLine(2);
    TPZVec <long> TopolPoint(1);
    
    //indice dos nos
    long id = 0;
    REAL valx;
    for(int xi = 0; xi < Qnodes/2; xi++)
    {
        valx = xi*Lx;
        Node[id].SetNodeId(id);
        Node[id].SetCoord(0 ,valx );//coord X
        Node[id].SetCoord(1 ,0. );//coord Y
        gmesh->NodeVec()[id] = Node[id];
        id++;
    }
    
    for(int xi = 0; xi < Qnodes/2; xi++)
    {
        valx = Lx - xi*Lx;
        Node[id].SetNodeId(id);
        Node[id].SetCoord(0 ,valx );//coord X
        Node[id].SetCoord(1 ,Ly);//coord Y
        gmesh->NodeVec()[id] = Node[id];
        id++;
    }
    
    //indice dos elementos
    id = 0;
    
//    if(triang_elements==true)
//    {
//        TopolTriang[0] = 0;
//        TopolTriang[1] = 1;
//        TopolTriang[2] = 3;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (id,TopolTriang,MatId,*gmesh);
//        id++;
//
//        TopolTriang[0] = 2;
//        TopolTriang[1] = 1;
//        TopolTriang[2] = 3;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,MatId,*gmesh);
//        id++;
//
//        TopolLine[0] = 0;
//        TopolLine[1] = 1;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC0,*gmesh);
//        id++;
//
//        TopolLine[0] = 2;
//        TopolLine[1] = 1;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC1,*gmesh);
//        id++;
//
//        TopolLine[0] = 3;
//        TopolLine[1] = 2;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC2,*gmesh);
//        id++;
//
//        TopolLine[0] = 3;
//        TopolLine[1] = 0;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC3,*gmesh);
//    }
//    else{
//        TopolQuad[0] = 0;
//        TopolQuad[1] = 1;
//        TopolQuad[2] = 2;
//        TopolQuad[3] = 3;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,MatId,*gmesh);
//        id++;
//
//        TopolLine[0] = 0;
//        TopolLine[1] = 1;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC0,*gmesh);
//        id++;
//
//        TopolLine[0] = 1;
//        TopolLine[1] = 2;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC1,*gmesh);
//        id++;
//
//        TopolLine[0] = 2;
//        TopolLine[1] = 3;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC2,*gmesh);
//        id++;
//
//        TopolLine[0] = 3;
//        TopolLine[1] = 0;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC3,*gmesh);
//    }
//
    gmesh->BuildConnectivity();
    
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout<<"\n\n Malha Geometrica Inicial\n ";
        gmesh->Print(sout);
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    
    return gmesh;
}


TPZGeoMesh *GMesh3(bool triang_elements){
    
    int Qnodes = 6;
    
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    gmesh->SetMaxNodeId(Qnodes-1);
    gmesh->NodeVec().Resize(Qnodes);
    TPZVec<TPZGeoNode> Node(Qnodes);
    
    TPZVec <long> TopolQuad(4);
    TPZVec <long> TopolTriang(3);
    TPZVec <long> TopolLine(2);
    //    TPZVec <long> TopolPoint(1);
    
    //indice dos nos
    long id = 0;
    REAL valx;
    REAL inix = -1.;
    REAL valy = -1.;
    for(int ix = 0; ix < Qnodes/2; ix++)
    {
        valx = inix + ix;
        Node[id].SetNodeId(id);
        Node[id].SetCoord(0 ,valx);//coord X
        Node[id].SetCoord(1 ,valy);//coord Y
        gmesh->NodeVec()[id] = Node[id];
        id++;
    }
    
    inix = 1.;
    valy = 1.;
    for(int ix = 0; ix < Qnodes/2; ix++)
    {
        valx = inix - ix;
        Node[id].SetNodeId(id);
        Node[id].SetCoord(0 ,valx);//coord X
        Node[id].SetCoord(1 ,valy);//coord Y
        gmesh->NodeVec()[id] = Node[id];
        id++;
    }
    
    //indice dos elementos
    id = 0;
    
//    if(triang_elements==true)
//    {
//        TopolTriang[0] = 0;
//        TopolTriang[1] = 1;
//        TopolTriang[2] = 5;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (id,TopolTriang,MatId,*gmesh);
//        id++;
//
//        TopolTriang[0] = 1;
//        TopolTriang[1] = 4;
//        TopolTriang[2] = 5;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,MatId,*gmesh);
//        id++;
//
//        TopolTriang[0] = 1;
//        TopolTriang[1] = 2;
//        TopolTriang[2] = 4;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,MatId,*gmesh);
//        id++;
//
//        TopolTriang[0] = 2;
//        TopolTriang[1] = 3;
//        TopolTriang[2] = 4;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,MatId,*gmesh);
//        id++;
//
//        TopolLine[0] = 0;
//        TopolLine[1] = 1;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC0,*gmesh);
//        id++;
//
//        TopolLine[0] = 1;
//        TopolLine[1] = 2;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC1,*gmesh);
//        id++;
//
//        TopolLine[0] = 2;
//        TopolLine[1] = 3;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC2,*gmesh);
//        id++;
//
//        TopolLine[0] = 3;
//        TopolLine[1] = 4;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC3,*gmesh);
//        id++;
//
//        TopolLine[0] = 4;
//        TopolLine[1] = 5;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC4,*gmesh);
//        id++;
//
//        TopolLine[0] = 5;
//        TopolLine[1] = 0;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC5,*gmesh);
//        id++;
//    }
//    else{
//        TopolQuad[0] = 0;
//        TopolQuad[1] = 1;
//        TopolQuad[2] = 4;
//        TopolQuad[3] = 5;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,MatId,*gmesh);
//        id++;
//
//        TopolQuad[0] = 1;
//        TopolQuad[1] = 2;
//        TopolQuad[2] = 3;
//        TopolQuad[3] = 4;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,MatId,*gmesh);
//        id++;
//
//        TopolLine[0] = 0;
//        TopolLine[1] = 1;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC0,*gmesh);
//        id++;
//
//        TopolLine[0] = 1;
//        TopolLine[1] = 2;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC1,*gmesh);
//        id++;
//
//        TopolLine[0] = 2;
//        TopolLine[1] = 3;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC2,*gmesh);
//        id++;
//
//        TopolLine[0] = 3;
//        TopolLine[1] = 4;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC3,*gmesh);
//        id++;
//
//        TopolLine[0] = 4;
//        TopolLine[1] = 5;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC4,*gmesh);
//        id++;
//
//        TopolLine[0] = 5;
//        TopolLine[1] = 0;
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC5,*gmesh);
//        id++;
//    }
    
    gmesh->BuildConnectivity();
    
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout<<"\n\n Malha Geometrica Inicial\n ";
        gmesh->Print(sout);
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    
    return gmesh;
}


void RefinamentoUnif(TPZGeoMesh* gmesh, int nDiv)
{
    for(int D = 0; D < nDiv; D++)
    {
        int nels = gmesh->NElements();
        for(int elem = 0; elem < nels; elem++)
        {
            TPZVec< TPZGeoEl * > filhos;
            TPZGeoEl * gel = gmesh->ElementVec()[elem];
            gel->Divide(filhos);
        }
    }
    //    gmesh->BuildConnectivity();
}

TPZCompMesh *CMeshFlux(ProblemConfig &config)
//TPZCompMesh *CMeshFlux(int pOrder,TPZGeoMesh *gmesh)
{
    /// criar materiais
    int dim = config.dimension;
    TPZGeoMesh *gmesh = config.gmesh;
    TPZMatPoisson3d *material;
    material = new TPZMatPoisson3d(MatId,dim);
    TPZMaterial * mat(material);
    material->NStateVariables();
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    
    cmesh->SetAllCreateFunctionsHDiv();
    
    
    cmesh->InsertMaterialObject(mat);
    cmesh->SetDefaultOrder(config.orderp);
    
    ///Inserir condicao de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    //bc materials
    TPZManVector<int,4> bcids(4,-1);
    for (auto matid : config.bcmaterialids) {
         TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 0.);
         int bctype = 0;
         TPZBndCond *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
         bc->TPZMaterial::SetForcingFunction(config.exact.Exact());
         cmesh->InsertMaterialObject(bc);
     }
    
//    if(IsHomogeneo==true)
//    {
//        ///Criar condicoes de contorno
//        TPZMaterial * BCond0 = material->CreateBC(mat, BC0,bcdirichlet, val1, val2);
//        TPZMaterial * BCond1 = material->CreateBC(mat, BC1,bcdirichlet, val1, val2);
//        TPZMaterial * BCond2 = material->CreateBC(mat, BC2,bcdirichlet, val1, val2);
//        TPZMaterial * BCond3 = material->CreateBC(mat, BC3,bcdirichlet, val1, val2);
//
////        TPZAutoPointer<TPZFunction<STATE> > bcmatNeumannAcima;
////        bcmatNeumannAcima = new TPZDummyFunction<STATE>(NeumannAcima);
////        BCond2->SetForcingFunction(bcmatNeumannAcima);
////
////        TPZAutoPointer<TPZFunction<STATE> > bcmatNeumannAbaixo;
////        bcmatNeumannAbaixo = new TPZDummyFunction<STATE>(NeumannAbaixo);
////        BCond0->SetForcingFunction(bcmatNeumannAbaixo);
////
////        TPZAutoPointer<TPZFunction<STATE> > bcmatDirichletEsquerda;
////        bcmatDirichletEsquerda = new TPZDummyFunction<STATE>(DirichletEsquerda);
////        BCond3->SetForcingFunction(bcmatDirichletEsquerda);
//
//        cmesh->InsertMaterialObject(BCond0);
//        cmesh->InsertMaterialObject(BCond1);
//        cmesh->InsertMaterialObject(BCond2);
//        cmesh->InsertMaterialObject(BCond3);
//    }else{
//        TPZMaterial * BCond0 = material->CreateBC(mat, BC0,bcdirichlet, val1, val2);
//        TPZMaterial * BCond1 = material->CreateBC(mat, BC1,bcdirichlet, val1, val2);
//        TPZMaterial * BCond2 = material->CreateBC(mat, BC2,bcdirichlet, val1, val2);
//        TPZMaterial * BCond3 = material->CreateBC(mat, BC3,bcdirichlet, val1, val2);
//        TPZMaterial * BCond4 = material->CreateBC(mat, BC4,bcdirichlet, val1, val2);
//        TPZMaterial * BCond5 = material->CreateBC(mat, BC5,bcdirichlet, val1, val2);
//
//        cmesh->InsertMaterialObject(BCond0);
//        cmesh->InsertMaterialObject(BCond1);
//        cmesh->InsertMaterialObject(BCond2);
//        cmesh->InsertMaterialObject(BCond3);
//        cmesh->InsertMaterialObject(BCond4);
//        cmesh->InsertMaterialObject(BCond5);
//    }
    
    cmesh->SetDimModel(config.dimension);
    cmesh->AutoBuild();//Ajuste da estrutura de dados computacional
    
    //#ifdef LOG4CXX
    //    if(logdata->isDebugEnabled())
    //    {
    //        std::stringstream sout;
    //        sout<<"\n\n Malha Computacional_1 Fluxo\n ";
    //        cmesh->Print(sout);
    //        LOGPZ_DEBUG(logdata,sout.str())
    //    }
    //#endif
    
    return cmesh;
}

TPZCompMesh *CMeshPressure(ProblemConfig &config)//(int pOrder,TPZGeoMesh *gmesh)
{
    TPZCompMesh *cmesh = new TPZCompMesh(config.gmesh);
    TPZMaterial *mat = 0;
    for (auto matid : config.materialids) {
        TPZMatPoisson3d *mix = new TPZMatPoisson3d(matid, cmesh->Dimension());
        if (!mat) mat = mix;
        cmesh->InsertMaterialObject(mix);
    }
    
    cmesh->SetDefaultOrder(config.orderp);
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    if(config.Iscontinuouspressure){
        cmesh->AutoBuild();
    }
    else{
        cmesh->ApproxSpace().CreateDisconnectedElements(true);
        cmesh->AutoBuild();
        int64_t n_connects = cmesh->NConnects();
        for (int64_t i = 0; i < n_connects; ++i) {
            cmesh->ConnectVec()[i].SetLagrangeMultiplier(1);
        }

    }
    
    
//    /// criar materiais
//    int dim = 2;
//    TPZMatPoisson3d *material;
//    material = new TPZMatPoisson3d(MatId,dim);
//    //material->NStateVariables();
//
//    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
//    cmesh->SetDimModel(dim);
//    TPZMaterial * mat(material);
//    cmesh->InsertMaterialObject(mat);
//    cmesh->SetDefaultOrder(pOrder);
//
//    ///Inserir condicao de contorno
//    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
//    if(IsHomogeneo==true)
//    {
//        TPZMaterial * BCond0 = material->CreateBC(mat, BC0,bcdirichlet, val1, val2);
//        TPZMaterial * BCond1 = material->CreateBC(mat, BC1,bcdirichlet, val1, val2);
//        TPZMaterial * BCond2 = material->CreateBC(mat, BC2,bcdirichlet, val1, val2);
//        TPZMaterial * BCond3 = material->CreateBC(mat, BC3,bcdirichlet, val1, val2);
//
//        cmesh->InsertMaterialObject(BCond0);
//        cmesh->InsertMaterialObject(BCond1);
//        cmesh->InsertMaterialObject(BCond2);
//        cmesh->InsertMaterialObject(BCond3);
//    }else{
//        TPZMaterial * BCond0 = material->CreateBC(mat, BC0,bcdirichlet, val1, val2);
//        TPZMaterial * BCond1 = material->CreateBC(mat, BC1,bcdirichlet, val1, val2);
//        TPZMaterial * BCond2 = material->CreateBC(mat, BC2,bcdirichlet, val1, val2);
//        TPZMaterial * BCond3 = material->CreateBC(mat, BC3,bcdirichlet, val1, val2);
//        TPZMaterial * BCond4 = material->CreateBC(mat, BC4,bcdirichlet, val1, val2);
//        TPZMaterial * BCond5 = material->CreateBC(mat, BC5,bcdirichlet, val1, val2);
//
//        cmesh->InsertMaterialObject(BCond0);
//        cmesh->InsertMaterialObject(BCond1);
//        cmesh->InsertMaterialObject(BCond2);
//        cmesh->InsertMaterialObject(BCond3);
//        cmesh->InsertMaterialObject(BCond4);
//        cmesh->InsertMaterialObject(BCond5);
//    }
//
//    if(IsContinuou==false)
//    {
//        cmesh->SetAllCreateFunctionsDiscontinuous();
//
//        //Ajuste da estrutura de dados computacional
//        cmesh->AutoBuild();
//
//        int ncon = cmesh->NConnects();
//        for(int i=0; i<ncon; i++)
//        {
//            TPZConnect &newnod = cmesh->ConnectVec()[i];
//            //newnod.SetPressure(true);
//            newnod.SetLagrangeMultiplier(1);
//        }
//
//        int nel = cmesh->NElements();
//        for(int i=0; i<nel; i++){
//            TPZCompEl *cel = cmesh->ElementVec()[i];
//            TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
//            celdisc->SetConstC(1.);
//            celdisc->SetCenterPoint(0, 0.);
//            celdisc->SetCenterPoint(1, 0.);
//            celdisc->SetCenterPoint(2, 0.);
//            celdisc->SetTrueUseQsiEta();
//            if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
//            {
//                if(fTriang==true) celdisc->SetTotalOrderShape();
//                else celdisc->SetTensorialShape();
//            }
//
//        }
//
//
//#ifdef PZDEBUG
//        int ncel = cmesh->NElements();
//        for(int i =0; i<ncel; i++){
//            TPZCompEl * compEl = cmesh->ElementVec()[i];
//            if(!compEl) continue;
//            TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
//            if(facel)DebugStop();
//
//        }
//#endif
//    }
//    else{
//        cmesh->SetAllCreateFunctionsContinuous();
//        //Ajuste da estrutura de dados computacional
//        cmesh->AutoBuild();
//    }
//
//    //#ifdef LOG4CXX
//    //    if(logdata->isDebugEnabled())
//    //    {
//    //        std::stringstream sout;
//    //        sout<<"\n\n Malha Computacional_2 pressure\n ";
//    //        cmesh->Print(sout);
//    //        LOGPZ_DEBUG(logdata,sout.str());
//    //    }
//    //#endif
    return cmesh;
}

TPZCompMesh *CMeshMixed(TPZVec<TPZCompMesh *> meshvec,ProblemConfig &config){//(TPZVec<TPZCompMesh *> meshvec,TPZGeoMesh * gmesh){
    


    TPZGeoMesh *gmesh= config.gmesh;
    //Creating computational mesh for multiphysic elements
    gmesh->ResetReference();
    TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    
    //criando material
    int dim = gmesh->Dimension();
    TPZMixedPoisson *material = new TPZMixedPoisson(MatId,dim);
    

    material->SetForcingFunction(config.exact.ForcingFunction());
    material->SetForcingFunctionExact(config.exact.Exact());
    
    //incluindo os dados do problema
//    REAL coefk = 1.;
//    material->SetPermeability(coefk);
//    REAL coefvisc = 1.;
//    material->SetViscosity(coefvisc);
    
    //permeabilidade
    if(IsHomogeneo==true){
        TPZFMatrix<REAL> Ktensor(3,3,0.);
        TPZFMatrix<REAL> InvK(3,3,0.);

        TPZFMatrix<REAL> K(3,3,0),invK(3,3,0);
        Ktensor.Identity();
        InvK.Identity();
        material->SetPermeabilityTensor(Ktensor,InvK);
    }
    else{
//        TPZAutoPointer<TPZFunction<STATE> > tensorK;
//        tensorK = new TPZDummyFunction<STATE>(PermeabilityTensor);
//        material->SetPermeabilityFunction(tensorK);
        //prcisamos implementar
        DebugStop();
    }
    
    if(IsStab==true){
        material->SetStabilizedMethod();
        material->SetStabilizationCoeficients(Delta1,Delta2);
    }
    if(IsStab==true && Useh2==true){
        material->SetHdois();
    }
    
    
    

    
    //inserindo o material na malha computacional
    TPZMaterial *mat(material);
    mphysics->InsertMaterialObject(mat);
    mphysics->SetDimModel(dim);
    //Criando condicoes de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    for (auto matid : config.bcmaterialids) {
        TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 0.);
        int bctype = 0;
        TPZBndCond *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
        bc->TPZMaterial::SetForcingFunction(config.exact.Exact());
        mphysics->InsertMaterialObject(bc);
    }
    
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    
    //Fazendo auto build
    mphysics->AutoBuild();
    mphysics->AdjustBoundaryElements();
    mphysics->CleanUpUnconnectedNodes();
    
    // Creating multiphysic elements into mphysics computational mesh
    TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
    TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
    
    return mphysics;
}

#include "TPZParFrontStructMatrix.h"
void ResolverSistema(TPZAnalysis &an, TPZCompMesh *fCmesh, int numthreads, bool direct)
{
    
    if (direct)
    {
        //        TPZSkylineStructMatrix full(fCmesh); //caso simetrico
        //        full.SetNumThreads(numthreads);
        //        an.SetStructuralMatrix(full);
        //        TPZStepSolver<STATE> step;
        //        step.SetDirect(ELDLt); //caso simetrico
        //        //  step.SetDirect(ELU);
        //        an.SetSolver(step);
        //        an.Run();
        
        TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(fCmesh);
        strmat.SetDecomposeType(ELDLt);
        strmat.SetNumThreads(numthreads);
        an.SetStructuralMatrix(strmat);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt); //caso simetrico
        //    step.SetDirect(ELU);
        an.SetSolver(step);
        an.Run();
    }
    else{
        
        TPZSkylineStructMatrix full(fCmesh); //caso simetrico
        full.SetNumThreads(numthreads);
        an.SetStructuralMatrix(full);
        TPZStepSolver<STATE> step;
        TPZAutoPointer<TPZMatrix<STATE> > Inverse = ComputeInverse(fCmesh);
        TPZStepSolver<STATE> precond(Inverse);
        precond.SetMultiply();
        //        step.SetGMRES(10, 40, precond, 1.0e-10, 0);
        step.SetCG(10, precond, 1.0e-10, 0);
        an.SetSolver(step);
        an.Run();
    }
    
    //Saida de Dados: solucao e  grafico no VT
    //    ofstream file("Solutout");
    //    an.Solution().Print("solution", file);    //Solution visualization on Paraview (VTK)
}

TPZFMatrix<STATE> * ComputeInverse(TPZCompMesh * mphysics)
{
    int neq = mphysics->NEquations();
    TPZFMatrix<STATE> * PreInverse =  new TPZFMatrix<STATE> (neq,neq,0.0);
    TPZSkylineStructMatrix skyl(mphysics);
    std::set<int> matids; // to be computed
    matids.insert(MatId);
    matids.insert(BC0);
    matids.insert(BC1);
    matids.insert(BC2);
    matids.insert(BC3);
    skyl.SetMaterialIds(matids);
    TPZFMatrix<STATE> rhsfrac;
    TPZAutoPointer<TPZGuiInterface> gui = new TPZGuiInterface;
    TPZAutoPointer<TPZMatrix<STATE> > matfrac = skyl.CreateAssemble(rhsfrac, gui);
    TPZFMatrix<STATE> oldmat = *matfrac.operator->();
    
    
    //#ifdef LOG4CXX
    //    if(logger->isDebugEnabled())
    //    {
    //        std::stringstream sout;
    //        sout.precision(20);
    //        matfrac->Print("K = ",sout,EMathematicaInput);
    //        LOGPZ_DEBUG(logger,sout.str());
    //    }
    //#endif
    
    matfrac->Inverse( * PreInverse,ELDLt);
    
    
    //#ifdef LOG4CXX
    //        if(logger->isDebugEnabled())
    //        {
    //            std::stringstream sout;
    //            sout.precision(20);
    ////            TPZFMatrix<> Iden(neq,neq,0.);
    ////            oldmat.Multiply(*PreInverse, Iden);
    //            PreInverse->Print("Kinv = ",sout,EMathematicaInput);
    ////            Iden.Print("Iden = ",sout,EMathematicaInput);
    //            LOGPZ_DEBUG(logger,sout.str());
    //        }
    //#endif
    
    
    return PreInverse;
    
}

void PosProcessMultph(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile){
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
    TPZManVector<std::string,10> scalnames(3), vecnames(4);
    vecnames[0]  = "Flux";
    vecnames[1]  = "GradFluxX";
    vecnames[2]  = "GradFluxY";
    vecnames[3]  = "ExactFlux";
    scalnames[0] = "Pressure";
    scalnames[1] = "DivFlux";
    scalnames[2] = "ExactPressure";
    
    
    const int dim = 2;
    int div =1;
    an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    an.PostProcess(div,dim);
    //    std::ofstream out("malha.txt");
    //    an.Print("nothing",out);
    
}

void SolExataMista(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux){
    
    solp.Resize(1, 0.);
    flux.Resize(3, 1.);
    flux(0,0)=flux(1,0)=flux(2,0)=0.;
    double x = pt[0];
    double y = pt[1];
    
    solp[0] = (-2./pi)*cos(pi*x)*exp(y/2.);
    flux(0,0)= -2*sin(pi*x)*exp(y/2.);
    flux(1,0)= (1./pi)*cos(pi*x)*exp(y/2.);
    flux(2,0)= (1.- 4.*pi*pi)/(2.*pi)*(cos(pi*x)*exp(y/2.));
}

void SolExataPressao(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux){
    
    solp.Resize(1, 0.);
    flux.Resize(3, 1.);
    flux(0,0)=flux(1,0)=flux(2,0)=0.;
    double x = pt[0];
    double y = pt[1];
    
    solp[0] = (-2./pi)*cos(pi*x)*exp(y/2.);
    flux(0,0)= 2*sin(pi*x)*exp(y/2.);
    flux(1,0)= -(1./pi)*cos(pi*x)*exp(y/2.);
    
}

void ForcingMista(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    double x = pt[0];
    double y = pt[1];
    disp[0]= (1.- 4.*pi*pi)/(2.*pi)*(cos(pi*x)*exp(y/2.));
}

void NeumannAbaixo(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
    double x = pt[0];
    //double y = pt[1];
    disp[0] = -cos(pi*x)/pi;
}

void NeumannAcima(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
    double x = pt[0];
    //double y = pt[1];
    disp[0] = 0.4087179842429694*cos(pi*x);
}

void DirichletEsquerda(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
    //    TPZFMatrix<STATE> flux;
    //    SolExataMista(pt, disp, flux);
    double y=pt[1];
    disp[0]=-(2./pi)*exp(y/2);
}

//meio heterogeneo
void SolFluxoHeter(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux){
    
    solp.Resize(1, 0.);
    flux.Resize(3, 1.);
    flux(0,0)=flux(1,0)=flux(2,0)=0.;
    REAL x = pt[0];
    REAL y = pt[1];
    
    if(x < 0.)
    {
        solp[0] = sin(y) + x*cos(y) + 2.*x*sin(y);
        flux(0,0)= -cos(y) - 2.*sin(y);//-KGradP[[1]]
        flux(1,0)= -cos(y) - 2.*x*cos(y) + x*sin(y);//-KGradP[[2]]
        flux(2,0)= x*cos(y) + sin(y) + 2.*x*sin(y);//div[-KgradP]
    }else
    {
        solp[0] = exp(x)*sin(y);
        flux(0,0) = -exp(x)*cos(y) - 2.*exp(x)*sin(y);//-KGradP[[1]]
        flux(1,0) = -2.*exp(x)*cos(y) - exp(x)*sin(y);//-KGradP[[2]]
        flux(2,0)= -2.*exp(x)*cos(y);//div[-KgradP]
    }
}

void SolPressaoHeter(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux){
    
    solp.Resize(1, 0.);
    flux.Resize(3, 1.);
    flux(0,0)=flux(1,0)=flux(2,0)=0.;
    REAL x = pt[0];
    REAL y = pt[1];
    
    if(x < 0.)
    {
        solp[0] = sin(y) + x*cos(y) + 2.*x*sin(y);
        flux(0,0)= cos(y) + 2.*sin(y);//dPdx
        flux(1,0)= cos(y) - x*sin(y) + 2.*x*cos(y);//dPdy
    }else
    {
        solp[0] = exp(x)*sin(y);
        flux(0,0) = exp(x)*sin(y);//dPdx
        flux(1,0) = exp(x)*cos(y);//dPdy
    }
    
}

void ForcingHeter(const TPZVec<REAL> &pt, TPZVec<STATE> &disp)
{
    REAL x = pt[0];
    REAL y = pt[1];
    
    if(x < 0.)
    {
        disp[0] = x*cos(y) + sin(y) + 2.*x*sin(y);
    }else
    {
        disp[0] = -2.*exp(x)*cos(y);
    }
}

void NeumannAbaixoXMenorZero(const TPZVec<REAL> &pt, TPZVec<STATE> &disp)
{
    REAL x = pt[0];
    if(x>0.) DebugStop();
    
    REAL normal[2] = {0.,-1.};
    TPZManVector<STATE> p(1);
    TPZFNMatrix<10,STATE> fluxo(2,1);
    SolFluxoHeter(pt,p,fluxo);
    
    disp.Resize(1);
    disp[0] = 0.5403023058681398 + 1.922075596544176*x;//fluxo(0,0)*normal[0]+fluxo(1,0)*normal[1];
}


void NeumannAbaixoXMaiorZero(const TPZVec<REAL> &pt, TPZVec<STATE> &disp)
{
    REAL x = pt[0];
    if(x<0.) DebugStop();
    
    REAL normal[2] = {0.,-1.};
    TPZManVector<STATE> p(1);
    TPZFNMatrix<10,STATE> fluxo(2,1);
    SolFluxoHeter(pt,p,fluxo);
    
    disp.Resize(1);
    disp[0] = 2.*exp(x)*cos(1.) - exp(x)*sin(1.);//fluxo(0,0)*normal[0]+fluxo(1,0)*normal[1];
    
}


void DirichletXIgualUm(const TPZVec<REAL> &pt, TPZVec<STATE> &disp)
{
    REAL x = pt[0];
    REAL y = pt[1];
    if(x != 1.) DebugStop();
    
    TPZManVector<STATE> p(1);
    TPZFNMatrix<10,STATE> fluxo(2,1);
    SolFluxoHeter(pt,p,fluxo);
    
    disp.Resize(1);
    disp[0] = 2.718281828459045*sin(y);//p[0];
}

void NeumannAcimaXMaiorZero(const TPZVec<REAL> &pt, TPZVec<STATE> &disp)
{
    REAL x = pt[0];
    if(x<0.) DebugStop();
    
    REAL normal[2] = {0.,1.};
    TPZManVector<STATE> p(1);
    TPZFNMatrix<10,STATE> fluxo(2,1);
    SolFluxoHeter(pt,p,fluxo);
    
    disp.Resize(1);
    disp[0] = -2.*exp(x)*cos(1.) - exp(x)*sin(1.);//fluxo(0,0)*normal[0]+fluxo(1,0)*normal[1];
}

void NeumannAcimaXMenorZero(const TPZVec<REAL> &pt, TPZVec<STATE> &disp)
{
    REAL x = pt[0];
    if(x>0.) DebugStop();
    
    REAL normal[2] = {0.,1.};
    TPZManVector<STATE> p(1);
    TPZFNMatrix<10,STATE> fluxo(2,1);
    SolFluxoHeter(pt,p,fluxo);
    
    disp.Resize(1);
    disp[0] = -0.23913362692838303*x - cos(1.);//fluxo(0,0)*normal[0]+fluxo(1,0)*normal[1];
}

void DirichletXIgualMenosUm(const TPZVec<REAL> &pt, TPZVec<STATE> &disp)
{
    REAL x = pt[0];
    REAL y = pt[1];
    if(x != -1.) DebugStop();
    
    TPZManVector<STATE> p(1);
    TPZFNMatrix<10,STATE> fluxo(2,1);
    SolFluxoHeter(pt,p,fluxo);
    
    disp.Resize(1);
    disp[0] = sin(y) - cos(y) - 2.*sin(y);//p[0];
}


void PermeabilityTensor(const TPZVec<REAL> &pt, TPZVec<STATE> &kabs, TPZFMatrix<STATE> &tensorK){
    
    
    if(IsHomogeneo){
        tensorK.Resize(4,2);
        
        tensorK(0,0)=1.;
        tensorK(0,1)=0.;
        tensorK(1,0)=0.;
        tensorK(1,1)=1.;
        
        //kinv
        tensorK(2,0)=1.;
        tensorK(2,1)=0.;
        tensorK(3,0)=0.;
        tensorK(3,1)=1.;
        
    }
    
    else{
        REAL x = pt[0];
        tensorK.Resize(4,2);
        kabs.Resize(1, 0.);
        
        if(x<0){
            //K
            tensorK(0,0)=1.;
            tensorK(0,1)=0.;
            tensorK(1,0)=0.;
            tensorK(1,1)=1.;
            
            //kinv
            tensorK(2,0)=1.;
            tensorK(2,1)=0.;
            tensorK(3,0)=0.;
            tensorK(3,1)=1.;
        }
        else{
            //K
            tensorK(0,0)=2.;
            tensorK(0,1)=1.;
            tensorK(1,0)=1.;
            tensorK(1,1)=2.;
            
            //kinv
            tensorK(2,0)=2./3.;
            tensorK(2,1)=-1./3.;
            tensorK(3,0)=-1./3.;
            tensorK(3,1)=2./3.;
        }
    }
}

void ErrorHDiv2(TPZCompMesh *hdivmesh, std::ostream &out)
{
    long nel = hdivmesh->NElements();
    int dim = hdivmesh->Dimension();
    TPZManVector<STATE,10> globerrors(10,0.);
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = hdivmesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel || gel->Dimension() != dim) {
            continue;
        }
        TPZManVector<REAL,10> elerror(10,0.);
        if(IsHomogeneo==true){
            cel->EvaluateError(SolExataMista, elerror, NULL);
        }else{
            cel->EvaluateError(SolFluxoHeter, elerror, NULL);
        }
        int nerr = elerror.size();
        for (int i=0; i<nerr; i++) {
            globerrors[i] += elerror[i]*elerror[i];
        }
        
    }
    // out << "Errors associated with HDiv space\n";
    out << "L2 Norm for flux = "    << sqrt(globerrors[1]) << std::endl;
    out << "L2 Norm for divergence = "    << sqrt(globerrors[2]) << std::endl;
    out << "Hdiv Norm for flux = "    << sqrt(globerrors[3])  <<std::endl;
    
}

void ErrorL22(TPZCompMesh *l2mesh, std::ostream &out,ProblemConfig &config)
{
//    TPZAnalysis an(l2mesh,false);
//    TPZManVector<REAL> errors(4,0.);
//    if(config.exact.Exact())
//    {
//
//        an.SetThreadsForError(0);
//        an.SetExact(config.exact.ExactSolution());
//        an.PostProcessError(errors,true);
//
//    }
        
    
    long nel = l2mesh->NElements();
    int dim = l2mesh->Dimension();
    TPZManVector<STATE,10> globerrors(10,0.);
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = l2mesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel || gel->Dimension() != dim) {
            continue;
        }
        TPZManVector<REAL,10> elerror(10,0.);
        cel->EvaluateError(config.exact.ExactSolution(), elerror, NULL);
//        if(IsHomogeneo==true){
//            cel->EvaluateError(SolExataPressao, elerror, NULL);
//        }else{
//            cel->EvaluateError(SolPressaoHeter, elerror, NULL);
//        }
        int nerr = elerror.size();
        globerrors.resize(nerr);
#ifdef LOG4CXX
        if (logger->isDebugEnabled()) {
            std::stringstream sout;
            sout << "L2 Error sq of element " << el << elerror[0]*elerror[0];
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        for (int i=0; i<nerr; i++) {
            globerrors[i] += elerror[i]*elerror[i];
        }

    }
    out << "\n";
    out << "Errors associated with L2 or H1 space\n";
    out << "H1 Norm = "    << sqrt(globerrors[0]) << "\n";
    out << "L2 Norm = "    << sqrt(globerrors[1]) << "\n";
    out << "Semi H1 Norm = " << sqrt(globerrors[2])<<"\n";
}

void ComputeFluxError(TPZCompMesh *cmesh, std::ostream &out){
    
    TPZManVector<REAL,2> errors(3,0.);
    
    int dimmesh = cmesh->Dimension();
    int nel = cmesh->NElements();
    int iel;
    for(iel=0; iel<nel; iel++)
    {
        TPZCompEl *cel = cmesh->ElementVec()[iel];
        if(!cel) continue;
        
        int dimcel = cel->Dimension();
        if(dimcel != dimmesh) continue;
        
        TPZGeoEl *gel = cel->Reference();
        TPZAutoPointer<TPZIntPoints> intrule = gel->CreateSideIntegrationRule(gel->NSides()-1, 2);
        TPZManVector<int,3> prevorder(dimcel), maxorder(dimcel,intrule->GetMaxOrder());
        intrule->GetOrder(prevorder);
        intrule->SetOrder(maxorder);
        
        TPZInterpolationSpace *sp = dynamic_cast<TPZInterpolationSpace *>(cel);
        TPZMaterialData data;
        sp->InitMaterialData(data);
        const int npoints = intrule->NPoints();
        TPZManVector<REAL,3> qsi(dimcel), xVec(3);
        
        for(int ip = 0; ip < npoints; ip++)
        {
            REAL weight;
            intrule->Point(ip,qsi,weight);
            sp->ComputeShape(qsi, data.x, data.jacobian, data.axes, data.detjac, data.jacinv, data.phi, data.dphi, data.dphix);
            weight *= fabs(data.detjac);
            sp->ComputeSolution(qsi,data);
            
            TPZManVector<STATE,2> flux(2,0.);
            STATE divfluxo;
            flux[0]=data.sol[0][0];
            flux[1]=data.sol[0][1];
            divfluxo =  data.dsol[0](0,0)+data.dsol[0](1,1);
            
            TPZManVector<STATE> uExato(1);
            TPZFNMatrix<100,STATE> duExato(2,1);
            gel->X(qsi,xVec);
            SolExataMista(xVec, uExato, duExato);
            
            
            TPZManVector<STATE,2> diff(2,0.);
            diff[0] = flux[0]- duExato(0,0);
            diff[1] = flux[1]- duExato(1,0);
            
            //erro L2 do fluxo
            errors[0] += weight*(diff[0]*diff[0] + diff[1]*diff[1]);
            
            //erro L2 do divergente do fluxo
            STATE diffDiv = abs(divfluxo - duExato(2,0));
            errors[1] += weight*diffDiv*diffDiv;
            
            //erro Hdiv para o fluxo
            //            errors[2] += errors[0] + errors[1];
        }
        intrule->SetOrder(prevorder);
    }
    
    //erro Hdiv para o fluxo
    errors[2] = errors[0] + errors[1];
    
    errors[0] = sqrt(errors[0]);
    errors[1] = sqrt(errors[1]);
    errors[2] = sqrt(errors[2]);
    
    out << "\n";
    out << "Erros associados ao fluxo na norma L2\n";
    out << "Norma L2  para fluxo = " << errors[0] << std::endl;
    out << "Norma L2 para divergente = " << errors[1] << std::endl;
    out << "Norma Hdiv para o fluxo = " << errors[2] << std::endl;
    
}///method

void ComputePressureError(TPZCompMesh *cmesh, std::ostream &out){
    
    TPZManVector<REAL,2> errors(3,0.);
    errors.Fill(0.);
    
    int dimmesh = cmesh->Dimension();
    int nel = cmesh->NElements();
    int iel;
    for(iel=0; iel<nel; iel++)
    {
        TPZCompEl *cel = cmesh->ElementVec()[iel];
        if(!cel) continue;
        
        int dimcel = cel->Dimension();
        if(dimcel != dimmesh) continue;
        
        TPZGeoEl *gel = cel->Reference();
        TPZAutoPointer<TPZIntPoints> intrule = gel->CreateSideIntegrationRule(gel->NSides()-1, 2);
        TPZManVector<int,3> prevorder(dimcel), maxorder(dimcel,intrule->GetMaxOrder());
        intrule->GetOrder(prevorder);
        intrule->SetOrder(maxorder);
        
        TPZInterpolationSpace *sp = dynamic_cast<TPZInterpolationSpace *>(cel);
        TPZMaterialData data;
        sp->InitMaterialData(data);
        const int npoints = intrule->NPoints();
        TPZManVector<REAL,3> qsi(dimcel), xVec(3);
        
        for(int ip = 0; ip < npoints; ip++)
        {
            REAL weight;
            intrule->Point(ip,qsi,weight);
            sp->ComputeShape(qsi, data.x, data.jacobian, data.axes, data.detjac, data.jacinv, data.phi, data.dphi, data.dphix);
            weight *= fabs(data.detjac);
            sp->ComputeSolution(qsi,data);
            
            TPZManVector<REAL,2> gradP(2,0.);
            REAL diffP;
            
            TPZManVector<STATE> uExato(1);
            TPZFNMatrix<100,STATE> duExato(2,1);
            gel->X(qsi,xVec);
            SolExataPressao(xVec, uExato, duExato);
            
            //erro L2 da pressao
            REAL solP = data.sol[0][0];
            diffP = solP - uExato[0];
            errors[1] += weight*(diffP*diffP);
            
            //erro semi H1 da pressao
            TPZFNMatrix<3,REAL> dsoldx;
            TPZFMatrix<REAL> dsoldaxes(2,1);
            dsoldaxes(0,0) = data.dsol[0][0];
            dsoldaxes(1,0) = data.dsol[0][1];
            
            TPZAxesTools<REAL>::Axes2XYZ(dsoldaxes, dsoldx, data.axes);
            TPZManVector<REAL,2> diffGrad(2,0.);
            diffGrad[0] = dsoldx(0,0)-duExato(0,0);
            diffGrad[1] = dsoldx(1,0)-duExato(1,0);
            errors[2] += weight*(diffGrad[0]*diffGrad[0] + diffGrad[1]*diffGrad[1]);
            
            //erro H1 para a pressao
            //errors[0] += errors[1] + errors[2];
        }
        intrule->SetOrder(prevorder);
    }
    
    //erro H1 para a pressao
    errors[0] = errors[1] + errors[2];
    
    errors[0] = sqrt(errors[0]);
    errors[1] = sqrt(errors[1]);
    errors[2] = sqrt(errors[2]);
    
    out << "\n";
    out << "Erros associados a pressao nas normas L2 e H1\n";
    out << "Norma H1 para a pressao = " << errors[0] << std::endl;
    out << "Norma L2 para a pressao = " << errors[1] << std::endl;
    out << "Norma semi-H1 para a pressao = " << errors[2] << std::endl;
    
    
}///method

void NEquationsCondensed(TPZCompMesh *cmesh, long &neqglob,long &neqcond, bool ismisto){
    
    long ncon = cmesh->NConnects();
    neqglob = 0;
    neqcond = 0;
    for(int i = 0; i< ncon; i++){
        TPZConnect &co  = cmesh->ConnectVec()[i];
        //if(co.HasDependency()) continue;
        if(co.HasDependency()  || co.IsCondensed() || !co.NElConnected() || co.SequenceNumber() == -1) continue;
        
        int nelc = co.NElConnected();
        if (nelc==0) DebugStop();
        
        int dofsize = co.NShape()*co.NState();
        neqglob += dofsize;
        
        //equacoes condensaveis
        if (nelc == 1){
            neqcond += dofsize;
        }
    }
    
    if(ismisto && IsContinuou==false){
        int nel2D =0;
        for(int i=0; i<cmesh->NElements(); i++){
            TPZCompEl *cel = cmesh->ElementVec()[i];
            if(!cel) continue;
            if(cel->Reference()->Dimension() == cmesh->Dimension()){
                nel2D++;
            }
        }
        neqcond = neqcond - nel2D;
    }
}

TPZMultiphysicsCompMesh *CreateHDivMesh(ProblemConfig &problem) {
   
    TPZMultiphysicsCompMesh *cmesh = new TPZMultiphysicsCompMesh(problem.gmesh);
    TPZMaterial *mat = NULL;
    TPZFMatrix<REAL> K(3,3,0),invK(3,3,0);
    K.Identity();
    invK.Identity();
    

//    K.Print(std::cout);
//    invK.Print(std::cout);
    
    for (auto matid : problem.materialids) {
        TPZMixedPoisson *mix = new TPZMixedPoisson(matid, cmesh->Dimension());
        mix->SetForcingFunction(problem.exact.ForcingFunction());
        mix->SetForcingFunctionExact(problem.exact.Exact());
        mix->SetPermeabilityTensor(K, invK);
        
        if (!mat) mat = mix;
        
        cmesh->InsertMaterialObject(mix);
    }
    for (auto matid : problem.bcmaterialids) {
        TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 0.);
        int bctype = 0;
        TPZBndCond *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
        bc->TPZMaterial::SetForcingFunction(problem.exact.Exact());
        cmesh->InsertMaterialObject(bc);
    }
    cmesh->ApproxSpace().SetAllCreateFunctionsMultiphysicElem();
   // std::set<int> matid;
//    matid.insert(1);
//    matid.insert(-1);
    TPZManVector<int> active(2,1);
    TPZManVector<TPZCompMesh *> meshvector(2,0);
    
    meshvector[0] = CMeshFlux(problem);
    meshvector[1] = CMeshPressure(problem);
    cmesh->BuildMultiphysicsSpace(active, meshvector);
    cmesh->LoadReferences();
//    bool keepmatrix = false;
//    bool keeponelagrangian = true;
//    TPZCompMeshTools::CreatedCondensedElements(cmesh, keeponelagrangian, keepmatrix);
    
    return cmesh;
}

void SolveStabilizedProblem(TPZCompMesh *cmesh,const ProblemConfig &config)
{
    
    TPZAnalysis an(cmesh,false);
    
    
#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(cmesh);
    strmat.SetNumThreads(0);
    //        strmat.SetDecomposeType(ELDLt);
    an.SetStructuralMatrix(strmat);
#else
    TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(cmesh_HDiv);
    strmat.SetNumThreads(0);
    //        TPZSkylineStructMatrix strmat3(cmesh_HDiv);
    //        strmat3.SetNumThreads(8);
#endif
    
    TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
    direct->SetDirect(ELDLt);
    an.SetSolver(*direct);
    delete direct;
    direct = 0;
    an.Assemble();
    an.Solve();//resolve o problema misto ate aqui
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("Pressure");
    scalnames.Push("ExactPressure");
    vecnames.Push("Flux");
    vecnames.Push("ExactFlux");
    
    int dim = config.gmesh->Dimension();
    
    std::stringstream sout;
    
    sout << config.dir_name << "/"  "StabProblem_"<<config.problemname<<"OrderP"<< config.orderp<<"OrderQ"<<config.orderq<<"Nref_"<<config.ndivisions<<".vtk";
    
    an.DefineGraphMesh(dim, scalnames, vecnames, sout.str());
    int resolution=2;
    an.PostProcess(resolution,dim);
    
    if(config.exact.Exact())
    {
        TPZManVector<REAL> errors(4,0.);
        an.SetThreadsForError(0);
        an.SetExact(config.exact.ExactSolution());
        an.PostProcessError(errors,false);
        
        //Erro
        
        ofstream myfile;
        myfile.open("MixedError.txt", ios::app);
        myfile << "\n\n Error for Mixed formulation " ;
        myfile << "\n-------------------------------------------------- \n";
        myfile << "Ndiv = " << config.ndivisions << " Order k = " << config.porder <<"\n";
        myfile << "Energy norm = " << errors[0] << "\n";//norma energia
        myfile << "error norm L2 = " << errors[1] << "\n";//norma L2
        myfile << "Semi norm H1 = " << errors[2] << "\n";//norma L2
        myfile.close();
        
    }
}
