//
//  TPZMixedStabilizedHdiv.hpp
//  StabilizedHdiv
//
//  Created by Denise De Siqueira on 20/01/20.
//

#ifndef TPZMixedStabilizedHdiv_hpp
#define TPZMixedStabilizedHdiv_hpp

#include <stdio.h>




/**
 This material implement the stabilized method proposed by Maicon Correa
 **/


#include <stdio.h>
#include "mixedpoisson.h"


class TPZMixedStabilizedHdiv : public TPZMixedPoisson
{

public:
    
    TPZMixedStabilizedHdiv(int matid, int dim);
    
    TPZMixedStabilizedHdiv();
    
    TPZMixedStabilizedHdiv(const TPZMixedStabilizedHdiv &copy);
    
    TPZMixedStabilizedHdiv(const TPZMixedPoisson &copy);
    
    virtual ~TPZMixedStabilizedHdiv();
    
    TPZMixedStabilizedHdiv &operator=(const TPZMixedStabilizedHdiv &copy);
    
    virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec) override;
    
    virtual int NEvalErrors() override {return 5;}
    
    virtual void  Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    
    int VariableIndex(const std::string &name) override;
    
    int NSolutionVariables(int var) override;
    
    void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout) override;
    
    virtual void Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors) override;
};

#endif /* TPZMixedStabilizedHdiv_hpp */
