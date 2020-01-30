//
//  TPZMixedStabilizedHdiv.cpp
//  StabilizedHdiv
//
//  Created by Denise De Siqueira on 20/01/20.
//

#include "TPZMixedStabilizedHdiv.h"
#include "pzaxestools.h"
#include "TPZAnalyticSolution.h"
#include "pzbndcond.h"

TPZMixedStabilizedHdiv::TPZMixedStabilizedHdiv(int matid, int dim): TPZMixedPoisson(matid,dim)
{
    
}

TPZMixedStabilizedHdiv::TPZMixedStabilizedHdiv():TPZMixedPoisson()
{
    
}

TPZMixedStabilizedHdiv::TPZMixedStabilizedHdiv(const TPZMixedStabilizedHdiv &copy):TPZMixedPoisson(copy)
{
    
}

TPZMixedStabilizedHdiv::TPZMixedStabilizedHdiv(const TPZMixedPoisson &copy):TPZMixedPoisson(copy)
{
    
}


TPZMixedStabilizedHdiv::~TPZMixedStabilizedHdiv()
{
    
}

TPZMixedStabilizedHdiv &TPZMixedStabilizedHdiv::operator=(const TPZMixedStabilizedHdiv &copy)
{
    TPZMixedPoisson::operator=(copy);
    return *this;
}

void TPZMixedStabilizedHdiv::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
///*
// A((u,p),(v,q))= (f,(v,q))
// A((u,p),(v,q)) = (InvKu, v) - (div v,p) - (div u,q)+(div u, div v)-(K gradp,gradq)
//(f,(v,q)) = -2*(f,q)+(f,div v)
// */
//
    STATE force = ff;
    if(fForcingFunction) {
        TPZManVector<STATE> res(1);
        fForcingFunction->Execute(datavec[1].x,res);
        force = res[0];
    }
    
    TPZFNMatrix<9,REAL> PermTensor;
    TPZFNMatrix<9,REAL> InvPermTensor;
    
    GetPermeability(datavec[1].x, PermTensor, InvPermTensor);
    
    // Setting the phis
    TPZFMatrix<REAL> &phiQ = datavec[0].phi;
    TPZFMatrix<REAL> &phip = datavec[1].phi;
    TPZFMatrix<REAL> &dphiQ = datavec[0].dphix;
    TPZFMatrix<REAL> &dphiP = datavec[1].dphix;
    TPZFNMatrix<9,REAL> dphiPXY(3,dphiP.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP, dphiPXY, datavec[1].axes);
    

    
    int phrq, phrp;
    phrp = phip.Rows();
    phrq = datavec[0].fVecShapeIndex.NElements();
    
    int nactive = 0;
    for (int i=0; i<datavec.size(); i++) {
        if (datavec[i].fActiveApproxSpace) {
            nactive++;
        }
    }
#ifdef PZDEBUG
    {
        if(phrp+phrq != ek.Rows())
        {
            DebugStop();
        }
    }
#endif
    //Calculate the matrix contribution for flux. Matrix A
    for(int iq=0; iq<phrq; iq++)
    {
        //ef(iq, 0) += 0.;
        int ivecind = datavec[0].fVecShapeIndex[iq].first;
        int ishapeind = datavec[0].fVecShapeIndex[iq].second;
        TPZFNMatrix<3,REAL> ivec(3,1,0.);
        for(int id=0; id<3; id++){
            ivec(id,0) = datavec[0].fNormalVec(id,ivecind);
        }
        
        //Inserindo termo de estabilizacao no termo de fonte
        REAL divqi = 0.;
//        if(fIsStabilized)
//        {
            //calculando div(qi)
            TPZFNMatrix<3,REAL> axesvec(3,1,0.);
            datavec[0].axes.Multiply(ivec,axesvec);
            for(int iloc=0; iloc<fDim; iloc++)
            {
                divqi += axesvec(iloc,0)*dphiQ(iloc,ishapeind);
            }
            ef(iq, 0) += weight*(0.5*divqi*force);
        //}
        
        TPZFNMatrix<3,REAL> ivecZ(3,1,0.);
        TPZFNMatrix<3,REAL> jvecZ(3,1,0.);
        for (int jq=0; jq<phrq; jq++)
        {
            TPZFNMatrix<3,REAL> jvec(3,1,0.);
            int jvecind = datavec[0].fVecShapeIndex[jq].first;
            int jshapeind = datavec[0].fVecShapeIndex[jq].second;
            
            for(int id=0; id<3; id++){
                jvec(id,0) = datavec[0].fNormalVec(id,jvecind);
            }
            
            //dot product between Kinv[u]v
            jvecZ.Zero();
            for(int id=0; id<fDim; id++){
                for(int jd=0; jd<fDim; jd++){
                    jvecZ(id,0) += InvPermTensor(id,jd)*jvec(jd,0);
                }
            }
            //jvecZ.Print("mat1 = ");
            REAL prod1 = ivec(0,0)*jvecZ(0,0) + ivec(1,0)*jvecZ(1,0) + ivec(2,0)*jvecZ(2,0);
            ek(iq,jq) += weight*phiQ(ishapeind,0)*phiQ(jshapeind,0)*prod1;
            
            
            //Inserindo termos de estabilizacao na matriz do fluxo
            //if(fIsStabilized==true)
            //{
                //termos de delta1
                //dot product between uKinv[v]
                ivecZ.Zero();
                for(int id=0; id<fDim; id++){
                    for(int jd=0; jd<fDim; jd++){
                        ivecZ(id,0) += InvPermTensor(id,jd)*ivec(jd,0);
                    }
                }
                //ivecZ.Print("mat2 = ");
                REAL prod2 = ivecZ(0,0)*jvec(0,0) + ivecZ(1,0)*jvec(1,0) + ivecZ(2,0)*jvec(2,0);
                ek(iq,jq) += (-1.)*weight*0.5*phiQ(ishapeind,0)*phiQ(jshapeind,0)*prod2;
                
                
                //termos de delta2: //dot product between divQu.divQv
                REAL divqj = 0.;
                TPZFNMatrix<3,REAL> axesvec(3,1,0.);
                datavec[0].axes.Multiply(jvec,axesvec);
                //calculando div(qj)
                for(int jloc=0; jloc<fDim; jloc++)
                {
                    divqj += axesvec(jloc,0)*dphiQ(jloc,jshapeind);
                }
                ek(iq,jq) += weight*0.5*divqi*divqj;
            //}
            
        }
    }
    
    
    // Coupling terms between flux and pressure. Matrix B
    for(int iq=0; iq<phrq; iq++)
    {
        int ivecind = datavec[0].fVecShapeIndex[iq].first;
        int ishapeind = datavec[0].fVecShapeIndex[iq].second;
        
        TPZFNMatrix<3,REAL> ivec(3,1,0.);
        for(int id=0; id<3; id++){
            ivec(id,0) = datavec[0].fNormalVec(id,ivecind);
            //ivec(1,0) = datavec[0].fNormalVec(1,ivecind);
            //ivec(2,0) = datavec[0].fNormalVec(2,ivecind);
        }
        TPZFNMatrix<3,REAL> axesvec(3,1,0.);
        datavec[0].axes.Multiply(ivec,axesvec);
        
        REAL divwq = 0.;
        for(int iloc=0; iloc<fDim; iloc++)
        {
            divwq += axesvec(iloc,0)*dphiQ(iloc,ishapeind);
        }
        for (int jp=0; jp<phrp; jp++) {
            
            REAL fact = (-1.)*weight*phip(jp,0)*divwq;
            // Matrix B
            ek(iq, phrq+jp) += fact;
            
            // Matrix B^T
            ek(phrq+jp,iq) += fact;
            
            
            //Inserindo termo de estabilizacao: delta1
            //if(fIsStabilized==true)
            //{
                //produto gardPu.Qv
                REAL dotVGradP = 0.;
                
                for(int k =0; k<fDim; k++)
                {
                    dotVGradP += ivec(k,0)*phiQ(ishapeind,0)*dphiP(k,jp);
                }
                
                REAL integration = (-1.)*weight*0.5*dotVGradP;
                
                // Estabilizacao delta1 na Matrix B
                ek(iq, phrq+jp) += integration;
                
                // Estabilizacao delta1 na Matrix BˆT
                ek(phrq+jp,iq) += integration;
            //}
        }
    }
    
    //termo fonte referente a equacao da pressao
    for(int ip=0; ip<phrp; ip++){
        ef(phrq+ip,0) += (-1.)*weight*force*phip(ip,0);
    }
    
    //Contribution for estabilization delta1 for gradPu*gradPv. Matrix D
   // if(fIsStabilized==true)
    //{
        //produto KgradPu x KgradPv
        TPZFNMatrix<3,REAL> dphiPuZ(dphiP.Rows(),dphiP.Cols(),0.);
        PermTensor.Multiply(dphiPXY, dphiPuZ);
        
        for(int ip=0; ip<phrp; ip++)
        {
            for(int jp=0; jp<phrp; jp++)
            {
                for(int k =0; k<3; k++)
                {
                    ek(phrq+ip, phrq+jp) += (-1.)*weight*0.5*dphiPuZ(k,ip)*dphiPXY(k,jp);
                }
                
            }
        }
    //}
}

void TPZMixedStabilizedHdiv::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
    
#ifdef PZDEBUG
    if (bc.Type() > 2 ) {
        std::cout << " Erro.!! Neste material utiliza-se apenas condicoes de Neumann e Dirichlet\n";
        DebugStop();
    }
#endif
    
    int dim = Dimension();
    
    TPZFMatrix<REAL>  &phiQ = datavec[0].phi;
    int phrq = phiQ.Rows();

    REAL v2 = bc.Val2()(0,0);
    REAL v1 = bc.Val1()(0,0);
    if(bc.HasForcingFunction())
    {
        TPZManVector<STATE> res(3);
        TPZFNMatrix<9,STATE> gradu(dim,1);
        bc.ForcingFunction()->Execute(datavec[0].x,res,gradu);
        if(bc.Type() == 0)
        {
            v2 = res[0];
        }
        else if(bc.Type() == 1 || bc.Type() == 2)
        {
            TPZFNMatrix<9,REAL> PermTensor, InvPermTensor;
            GetPermeability(datavec[0].x, PermTensor, InvPermTensor);
            REAL normflux = 0.;
            for(int i=0; i<3; i++)
            {
                for(int j=0; j<dim; j++)
                {
                    normflux += datavec[0].normal[i]*PermTensor(i,j)*gradu(j,0);
                }
            }
            v2 = -normflux;
            if(bc.Type() ==2)
            {
                v2 = -res[0]+v2/v1;
            }
        }
        else
        {
            DebugStop();
        }
    }
    else
    {
        v2 = bc.Val2()(0,0);
    }

    switch (bc.Type()) {
        case 0 :        // Dirichlet condition
            //primeira equacao
            for(int iq=0; iq<phrq; iq++)
            {
                //the contribution of the Dirichlet boundary condition appears in the flow equation
                ef(iq,0) += (-1.)*v2*phiQ(iq,0)*weight;
            }
            break;
            
        case 1 :            // Neumann condition
            //primeira equacao
            for(int iq=0; iq<phrq; iq++)
            {
                ef(iq,0)+= gBigNumber*v2*phiQ(iq,0)*weight;
                for (int jq=0; jq<phrq; jq++) {
                    
                    ek(iq,jq)+= gBigNumber*phiQ(iq,0)*phiQ(jq,0)*weight;
                }
            }
            break;
        
        case 2 :            // mixed condition
            for(int iq = 0; iq < phrq; iq++) {
                
                ef(iq,0) += v2*phiQ(iq,0)*weight;
                for (int jq = 0; jq < phrq; jq++) {
                    ek(iq,jq) += weight/v1*phiQ(iq,0)*phiQ(jq,0);
                }
            }
            break;
    }
    
}


int TPZMixedStabilizedHdiv::VariableIndex(const std::string &name)
{
    if(name == "FluxFem") return 40;
    if(name == "FluxExact") return 41;
    if(name == "DivFluxFem") return 42;
    if(name == "DivFluxExact") return 43;
    if(name == "PressureFem") return 44;
    if(name == "PressureExact") return 45;
    if(name == "POrderPressure") return 46;
    if(name == "GradPressureFem") return 47;
    
    return -1;
}

int TPZMixedStabilizedHdiv::NSolutionVariables(int var)
{
    switch (var) {
        case 40:
        case 41:
            return 3;
            break;
        case 42:
        case 43:
        case 44:
        case 45:
        case 46:
            return 1;
            break;
        case 47:
            return 3;
            break;
        default:
            DebugStop();
            break;
    }
    return 0;
}

void TPZMixedStabilizedHdiv::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout)
{
    /**
     datavec[0] Hdiv mesh fem
     datavec[1] L2 mesh fem
     **/
    
    TPZFNMatrix<9,REAL> PermTensor = fTensorK;
    TPZFNMatrix<9,REAL> InvPermTensor = fInvK;
    
    
    TPZManVector<STATE,2> pressexact(1,0.);
    TPZFNMatrix<9,STATE> gradu(3,1,0.), fluxinv(3,1);
    
    TPZVec<STATE> divsigma_exact(1);
    if(fForcingFunctionExact)
    {
        this->fForcingFunctionExact->Execute(datavec[0].x, pressexact,gradu);
        this->fForcingFunction->Execute(datavec[0].x,divsigma_exact);
        
    }
    PermTensor.Multiply(gradu,fluxinv);
    
    //calculo do gradiente da pressao
    TPZFNMatrix<3,REAL> dsoldx;
    TPZFMatrix<REAL> dsoldaxes(fDim,1);
    for (int i=0; i<fDim; i++) dsoldaxes(i,0) = datavec[1].dsol[0][i];
    TPZAxesTools<REAL>::Axes2XYZ(dsoldaxes, dsoldx, datavec[1].axes);
    
    switch (var)
    {
        case 40://FluxFem
            for(int i=0; i<fDim; i++) Solout[i] = datavec[0].sol[0][i];
            break;
        case 41://FluxExact
            for(int i=0; i<fDim; i++) Solout[i] = -fluxinv(i);
            break;
        case 42://DivFluxFem
            Solout[0] = datavec[0].divsol[0][0];
            break;
        case 43://DivFluxExact
             Solout[0] = divsigma_exact[0];
            break;
        case 44://PressureFem
            Solout[0] = datavec[1].sol[0][0];
            break;
        case 45://PressureExact
            Solout[0] = pressexact[0];
            break;
        case 46://POrderPressure
            Solout[0] = datavec[1].p;
            break;
        case 47:
            Solout[0]=0.;
            for (int i=0; i<fDim; i++) Solout[i] = dsoldx(i,0);
            break;
        default:
            DebugStop();
    }
}

void TPZMixedStabilizedHdiv::Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors)
{
    /**
     datavec[0] Flux
     datavec[1] L2 mesh,
     
     error[0] - eror em norma L2 para a variavel pressao
     error[1] - erro em semi norma H1 para a variavel pressao
     error[2] - erro em norma H1
     error[3] - eror em norma L2 para a variavel fluxo
     error[4] - eror em norma L2 para a variavel divergente do fluxo
     **/
    
    errors.Resize(NEvalErrors());
    errors.Fill(0.0);
    
    TPZManVector<STATE,3> fluxfem(3), pressurefem(1);
    fluxfem = data[0].sol[0];
    pressurefem[0] = data[1].sol[0][0];
    STATE divsigmafem = data[0].divsol[0][0];
    
    TPZVec<STATE> divsigma_exact(1);
    if(this->fForcingFunctionExact){
        this->fForcingFunctionExact->Execute(data[0].x,u_exact,du_exact);
        this->fForcingFunction->Execute(data[0].x,divsigma_exact);
    }
    
    REAL errodivsigma = 0.;
    errodivsigma = (divsigma_exact[0] - divsigmafem)*(divsigma_exact[0] - divsigmafem);
    
    TPZFNMatrix<9,REAL> PermTensor = fTensorK;
    TPZFNMatrix<9,REAL> InvPermTensor = fInvK;
    
    
    TPZFNMatrix<3,REAL> fluxexactneg;
    
    //sigma=-K grad(u)
    
    {
        TPZFNMatrix<3,REAL> gradpressure_exact(3,1);
        for (int i=0; i<fDim; i++) {
            gradpressure_exact(i,0) = du_exact[i];
        }
        PermTensor.Multiply(gradpressure_exact,fluxexactneg);
    }
    
    
    REAL innerexact = 0.;
    //REAL innerestimate = 0.;
    for (int i=0; i<fDim; i++) {
        for (int j=0; j<fDim; j++) {
            innerexact += (fluxfem[i]+fluxexactneg(i,0))*InvPermTensor(i,j)*(fluxfem[j]+fluxexactneg(j,0));//Pq esta somando: o fluxo fem esta + e o exato -
        }
    }
    errors[0] = (pressurefem[0]-u_exact[0])*(pressurefem[0]-u_exact[0]);//exact error pressure
    
    TPZManVector<STATE,3> gradpressure(3,0.);
    this->Solution(data,VariableIndex("GradPressureFem"), gradpressure);
    REAL errogradpressure = 0.;
    for (int i=0; i<fDim; i++) {
//        for (int j=0; j<fDim; j++) {
//            errogradpressure += (gradpressure[i] - du_exact(i,0))*PermTensor(i,j)*(gradpressure[i] - du_exact(j,0));
//        }
        errogradpressure += (gradpressure[i] - du_exact(i,0))*(gradpressure[i] - du_exact(i,0));
    }

    errors[1] = errogradpressure;
    errors[2] = errors[0] + errors[1];
    
    errors[3] = innerexact;//error flux exact
    errors[4] = errodivsigma; //||div(sigma) - div(sigma_h)||
}
