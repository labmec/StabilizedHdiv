//
//  TPZMixedStabilizedHdiv.hpp
//  StabilizedHdiv
//
//  Created by Denise De Siqueira on 20/01/20.
//

#ifndef TPZMixedStabilizedHdiv_hpp
#define TPZMixedStabilizedHdiv_hpp

#include <stdio.h>

#endif /* TPZMixedStabilizedHdiv_hpp */


/**
 This material implement the stabilized method proposed by Maicon Correa
 **/


#include <stdio.h>
#include "mixedpoisson.h"


class TPZMixedStabilizedHdiv : public TPZMixedPoisson
{

public:
    
    
    TPZMixedStabilizedHdiv(int matid, int dim);
    //{
//        fdelta1 = 0.5;
//        fdelta2 = 0.5;
//        fh2 = 1.;
//    }
    
    TPZMixedStabilizedHdiv();
    
    TPZMixedStabilizedHdiv(const TPZMixedStabilizedHdiv &copy);
    
    TPZMixedStabilizedHdiv(const TPZMixedPoisson &copy);
    
    virtual ~TPZMixedStabilizedHdiv();
    
    TPZMixedStabilizedHdiv &operator=(const TPZMixedStabilizedHdiv &copy);
    //{
//        fdelta1 = copy.fdelta1;
//        fdelta2 = copy.fdelta2;
//        fh2 = copy.fh2;
//    }
    
    virtual void  Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;

};
