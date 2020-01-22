//
//  TPZMixedStabilizedHdiv.cpp
//  StabilizedHdiv
//
//  Created by Denise De Siqueira on 20/01/20.
//

#include "TPZMixedStabilizedHdiv.h"
#include "pzaxestools.h"

void TPZMixedStabilizedHdiv::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
/*
 A((u,p),(v,q))= (f,(v,q))
 A((u,p),(v,q)) = (InvKu, v) - (div v,p) - (div u,q)+(div u, div v)-(K gradp,gradq)
(f,(v,q)) = -2*(f,q)+(f,div v)
 */
    
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
    
    
   // REAL &faceSize = datavec[0].HSize;
    
  //      fh2 = faceSize*faceSize;
    
    
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
    if(nactive == 4)
    {
        int phrgb = datavec[2].phi.Rows();
        int phrub = datavec[3].phi.Rows();
        if(phrp+phrq+phrgb+phrub != ek.Rows())
        {
            DebugStop();
        }
    }else
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
        //if(fIsStabilized)
       // {
            //calculando div(qi)
            TPZFNMatrix<3,REAL> axesvec(3,1,0.);
            datavec[0].axes.Multiply(ivec,axesvec);
            for(int iloc=0; iloc<fDim; iloc++)
            {
                divqi += axesvec(iloc,0)*dphiQ(iloc,ishapeind);
            }
            ef(iq, 0) += weight*(divqi*force);
       // }
        
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
//            jvecZ.Zero();
//            for(int id=0; id<fDim; id++){
//                for(int jd=0; jd<fDim; jd++){
//                    jvecZ(id,0) += InvPermTensor(id,jd)*jvec(jd,0);
//                }
//            }
//            //jvecZ.Print("mat1 = ");
//            REAL prod1 = ivec(0,0)*jvecZ(0,0) + ivec(1,0)*jvecZ(1,0) + ivec(2,0)*jvecZ(2,0);
//            ek(iq,jq) += fvisc*weight*phiQ(ishapeind,0)*phiQ(jshapeind,0)*prod1;
            
            
            //Inserindo termos de estabilizacao na matriz do fluxo
           // if(fIsStabilized==true)
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
                ek(iq,jq) += (-1.)*weight*phiQ(ishapeind,0)*phiQ(jshapeind,0)*prod2;
                
                
                //termos de delta2: //dot product between divQu.divQv
                REAL divqj = 0.;
                TPZFNMatrix<3,REAL> axesvec(3,1,0.);
                datavec[0].axes.Multiply(jvec,axesvec);
                //calculando div(qj)
                for(int jloc=0; jloc<fDim; jloc++)
                {
                    divqj += axesvec(jloc,0)*dphiQ(jloc,jshapeind);
                }
                ek(iq,jq) += weight*divqi*divqj;
           // }
            
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
          //  if(fIsStabilized==true)
          //  {
                //produto gardPu.Qv
                REAL dotVGradP = 0.;
    
                for(int k =0; k<fDim; k++)
                {
                    dotVGradP += ivec(k,0)*phiQ(ishapeind,0)*dphiP(k,jp);
                }
                
                REAL integration = (-1.)*weight*dotVGradP;
                
                // Estabilizacao delta1 na Matrix B
                ek(iq, phrq+jp) += integration;
                
                // Estabilizacao delta1 na Matrix BˆT
                ek(phrq+jp,iq) += integration;
           // }
        }
    }
    
    //termo fonte referente a equacao da pressao
    for(int ip=0; ip<phrp; ip++){
        ef(phrq+ip,0) += (-2.)*weight*force*phip(ip,0);
    }
    
    //Contribution for estabilization delta1 for gradPu*gradPv. Matrix D
   // if(fIsStabilized==true)
   // {
        //produto KgradPu x KgradPv
        TPZFNMatrix<3,REAL> dphiPuZ(dphiP.Rows(),dphiP.Cols(),0.);
        PermTensor.Multiply(dphiPXY, dphiPuZ);
        
        for(int ip=0; ip<phrp; ip++)
        {
            for(int jp=0; jp<phrp; jp++)
            {
                for(int k =0; k<3; k++)
                {
                    ek(phrq+ip, phrq+jp) += (-1.)*weight*dphiPuZ(k,ip)*dphiPXY(k,jp);
                }
                
            }
        }
  //  }
    if(nactive == 4)
    {
        for(int ip=0; ip<phrp; ip++)
        {
            ek(phrq+ip,phrq+phrp) += phip(ip,0)*weight;
            ek(phrq+phrp,phrq+ip) += phip(ip,0)*weight;
        }
        ek(phrp+phrq+1,phrq+phrp) += -weight;
        ek(phrq+phrp,phrp+phrq+1) += -weight;
    }

    
}


