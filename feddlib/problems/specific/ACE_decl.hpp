#ifndef ACE_decl_hpp
#define ACE_decl_hpp
#include "feddlib/problems/abstract/Problem.hpp"
#include "feddlib/core/FE/ace_layer/ace_layer.hpp"

/*
The Ace -> Specific Problem.
It holds a pointer to the AceGenElement which is initialized by the constructor.
For the future, this class is supposed to handle a variable number of "DomainConstPtr_Type" (fields)
and also loses its restriction to a block matrix.
*/

namespace FEDD {
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class ACE : public Problem<SC,LO,GO,NO> {
    
public:
    
    
    typedef Problem<SC,LO,GO,NO> Problem_Type;
    typedef typename Problem_Type::Matrix_Type Matrix_Type;
    typedef typename Problem_Type::MatrixPtr_Type MatrixPtr_Type;
    typedef typename Problem_Type::Map_Type Map_Type;    
    typedef typename Problem_Type::MapPtr_Type MapPtr_Type;
    typedef typename Problem_Type::MapConstPtr_Type MapConstPtr_Type;
    
    typedef typename Problem_Type::BlockMatrix_Type BlockMatrix_Type;
    
    typedef typename Problem_Type::MultiVector_Type MultiVector_Type;
    typedef typename Problem_Type::MultiVectorPtr_Type MultiVectorPtr_Type;
    typedef typename Problem_Type::MultiVectorConstPtr_Type MultiVectorConstPtr_Type;

    typedef typename Problem_Type::DomainConstPtr_Type DomainConstPtr_Type;
    typedef typename Problem_Type::CommConstPtr_Type CommConstPtr_Type;
    
    ACE( const DomainConstPtr_Type &domainVelocity, std::string FETypeVelocity, const DomainConstPtr_Type &domainPressure, std::string FETypePressure, ParameterListPtr_Type parameterList, Teuchos::RCP<AceGenElement> AceElmt);
    ~ACE();
    
    virtual void info();
    
    virtual void assemble( std::string type ) const;
    
//    virtual void assembleExternal( std::string type );

    virtual void getValuesOfInterest( vec_dbl_Type& values ){};
    
    virtual void computeValuesOfInterestAndExport() {};
    
  private:
    mutable MultiVectorPtr_Type u_repNewton_;
    mutable MultiVectorPtr_Type p_repNewton_;
    mutable MultiVectorPtr_Type u_repTime_;
    mutable MultiVectorPtr_Type p_repTime_;
    /*####################*/

  public:
    Teuchos::RCP<AceGenElement> AceElmt;

};
}
#endif
