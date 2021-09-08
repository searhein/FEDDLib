#ifndef ACE_def_hpp
#define ACE_def_hpp
#include "ACE_decl.hpp"

namespace FEDD {
void ZeroDirichlet(double* x, double* res, double t, double* parameters){
    
    res[0] = 0.;
    
    return;
}
    
template<class SC,class LO,class GO,class NO>
ACE<SC,LO,GO,NO>::ACE(const DomainConstPtr_Type &domainVelocity, std::string FETypeVelocity, const DomainConstPtr_Type &domainPressure, std::string FETypePressure, ParameterListPtr_Type parameterList):
Problem<SC,LO,GO,NO>(parameterList, domainVelocity->getComm())
{

    this->addVariable( domainVelocity , FETypeVelocity , "u" , domainVelocity->getDimension());
    this->addVariable( domainPressure , FETypePressure , "p" , 1);
    this->dim_ = this->getDomain(0)->getDimension();
}

template<class SC,class LO,class GO,class NO>
ACE<SC,LO,GO,NO>::~ACE(){

}

template<class SC,class LO,class GO,class NO>
void ACE<SC,LO,GO,NO>::info(){
    this->infoProblem();
}
    
template<class SC,class LO,class GO,class NO>
void ACE<SC,LO,GO,NO>::assemble( std::string type ) const{
    
    if ( type == "FirstAssemble") {
        
        if (u_repNewton_.is_null()){
            u_repNewton_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated() )) ;
            u_repNewton_->putScalar(0.);
        }
        if (p_repNewton_.is_null()){
            p_repNewton_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(1)->getMapRepeated() )) ;
            p_repNewton_->putScalar(0.);
        }
        
        if (u_repTime_.is_null()){
            u_repTime_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated() )) ;
            u_repTime_->putScalar(0.);
        }
        if (p_repTime_.is_null()){
            p_repTime_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(1)->getMapRepeated() )) ;
            p_repTime_->putScalar(0.);
        }

    }
    
    if (type == "SetSolutionInTime") {
        TEUCHOS_TEST_FOR_EXCEPTION( u_repTime_.is_null(), std::runtime_error, "u_repTime_ not initialized.");
        TEUCHOS_TEST_FOR_EXCEPTION( p_repTime_.is_null(), std::runtime_error, "p_repTime_ not initialized.");
        MultiVectorConstPtr_Type u = this->solution_->getBlock(0);
        MultiVectorConstPtr_Type p = this->solution_->getBlock(1);
        u_repTime_->importFromVector(u, true);
        p_repTime_->importFromVector(p, true);
    }
    else if (type == "SetSolutionNewton") {
        TEUCHOS_TEST_FOR_EXCEPTION( u_repNewton_.is_null(), std::runtime_error, "u_repNewton_ not initialized.");
        TEUCHOS_TEST_FOR_EXCEPTION( p_repNewton_.is_null(), std::runtime_error, "p_repNewton_ not initialized.");

        MultiVectorConstPtr_Type u = this->solution_->getBlock(0);
        MultiVectorConstPtr_Type p = this->solution_->getBlock(1);
        u_repNewton_->importFromVector(u, true);
        p_repNewton_->importFromVector(p, true);
    }
    else if( type == "FirstAssemble" || type == "Assemble" || type == "OnlyUpdate" || type == "AssembleAndUpdate"){
        if (this->verbose_)
            std::cout << "-- External assembly ... " << std::flush;
    
        MapConstPtr_Type mapRepeatedConst1 = this->getDomain(0)->getMapRepeated();
        MapConstPtr_Type mapRepeatedConst2 = this->getDomain(1)->getMapRepeated();
        MapPtr_Type mapRepeated1 = Teuchos::rcp_const_cast<Map_Type>(mapRepeatedConst1);
        MapPtr_Type mapRepeated2 = Teuchos::rcp_const_cast<Map_Type>(mapRepeatedConst2);
                
        MatrixPtr_Type A00(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
        MatrixPtr_Type A01(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(1)->getDimension() * this->getDomain(1)->getApproxEntriesPerRow() ) );
        MatrixPtr_Type A10(new Matrix_Type( this->getDomain(1)->getMapUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
        MatrixPtr_Type A11(new Matrix_Type( this->getDomain(1)->getMapUnique(), this->getDomain(1)->getDimension() * this->getDomain(1)->getApproxEntriesPerRow() ) );
        
        MultiVectorPtr_Type F0(new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated(), 1 ) );
        MultiVectorPtr_Type F1(new MultiVector_Type( this->getDomain(1)->getMapRepeated(), 1 ) );
        
        bool update(type == "FirstAssemble" || type == "Assemble" || type == "AssembleAndUpdate");
        bool updateHistory(type == "OnlyUpdate" || type == "AssembleAndUpdate"); //for AceGen History variables

        this->feFactory_->assemblyAceGenACE(A00,
                                            A01,
                                            A10,
                                            A11,
                                            F0,
                                            F1,
                                            mapRepeated1,
                                            mapRepeated2,
                                            this->parameterList_,
                                            u_repNewton_,
                                            p_repNewton_,
                                            u_repTime_,
                                            p_repTime_,
                                            update,
                                            updateHistory);

        if (update) {
            this->system_->addBlock( A00, 0, 0 );
            this->system_->addBlock( A01, 0, 1 );
            this->system_->addBlock( A10, 1, 0 );
            this->system_->addBlock( A11, 1, 1 );

            MultiVectorPtr_Type F0Unique = Teuchos::rcp(new MultiVector_Type( this->getDomain(0)->getMapVecFieldUnique() ) );
            MultiVectorPtr_Type F1Unique = Teuchos::rcp(new MultiVector_Type( this->getDomain(1)->getMapUnique() ) );
            F0Unique->exportFromVector( F0, false, "Add" );
            F1Unique->exportFromVector( F1, false, "Add" );
            
            this->rhs_->addBlock( F0Unique, 0 );
            this->rhs_->addBlock( F1Unique, 1 );
            
        }
        
        
        if (this->verbose_)
            std::cout << "done -- " << std::endl;
    }
}

}
#endif
