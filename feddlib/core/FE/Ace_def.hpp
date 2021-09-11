/*
The new assembly routine that will call the AceGenElement
*/

#include "feddlib/core/FE/ace_layer/ace_layer.hpp"

template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::assemblyAceGenACE(    MatrixPtr_Type &A00,
                                            MatrixPtr_Type &A01,
                                            MatrixPtr_Type &A10,
                                            MatrixPtr_Type &A11,
                                            MultiVectorPtr_Type &F0,
                                            MultiVectorPtr_Type &F1,
                                            MapPtr_Type &mapRepeated1,
                                            MapPtr_Type &mapRepeated2,
                                            ParameterListPtr_Type parameterList,
                                            Teuchos::RCP<AceGenElement> AceElmt,
                                            MultiVectorPtr_Type u_repeatedNewton,
                                            MultiVectorPtr_Type p_repeatedNewton,
                                            MultiVectorPtr_Type u_repeatedTime,
                                            MultiVectorPtr_Type p_repeatedTime,
                                            bool update,
                                            bool updateHistory
                                            )
{
    
    // testen des elements: das hat funktioniert !!
    // int nonodes = AceElmt -> get_NoNodes();
    // int nodofs = AceElmt ->  get_NoDofs();
    // std::cout << "ELEMENTNODES: " << nonodes << std::endl;
    // std::cout << "ELEMENTDOFS : " << nodofs << std::endl;
    

    std::string tpmType = parameterList->sublist("Parameter").get("TPM Type","Biot");
    
    int dim = domainVec_[0]->getDimension();
    int idata = 1; //= we should init this
    int ic = -1; int ng = -1;
    
    //ed.hp:history previous (timestep); previous solution (velocity and acceleration)
    //ed.ht:? same length as hp
    ElementsPtr_Type elements1 = domainVec_[0]->getElementsC();
    ElementsPtr_Type elements2 = domainVec_[1]->getElementsC();
    
    int sizeED = 24; /* 2D case for P2 elements:
                      12 velocities, 12 accelerations (2 dof per P2 node)
                      */
    if (dim==3)
        sizeED = 60;/* 3D case for P2 elements:
                       30 velocities, 30 accelerations (3 dof per P2 node)
                    */
    // there is a element data field already, ed_ which holds the history data
    // see decl
    // this way the data are cept through multiple calls of this assembly
    if (ed_.size()==0){
        for (UN T=0; T<elements1->numberElements(); T++)
            ed_.push_back( Teuchos::rcp(new DataElement( sizeED )) );
    }
    
    ///---------------------------------------------------------------------------------------------------------------------------------------
    //              erstellen eines vectors von element spec ... das sind nur temporaere daten 
    //              er liest aber spaeter tatsaechlich die richtigen daten ein
    ///---------------------------------------------------------------------------------------------------------------------------------------
    std::vector<ElementSpec> es_vec( parameterList->sublist("Parameter").get("Number of materials",1) , ElementSpec());
    

    ///---------------------------------------------------------------------------------------------------------------------------------------
    //              Einlesen von "material parameter"
    ///---------------------------------------------------------------------------------------------------------------------------------------
    
    // erstellen eines parameter vectors fuer jedes domain
    vec2D_dbl_Type dataVec( parameterList->sublist("Parameter").get("Number of materials",1), vec_dbl_Type(6,0.) );

    for (int i=0; i<dataVec.size(); i++) {
        if (tpmType == "Biot") {
            if (dim==2) {
                dataVec[i][0] = parameterList->sublist("Timestepping Parameter").get("Newmark gamma",0.5);
                dataVec[i][1] = parameterList->sublist("Timestepping Parameter").get("Newmark beta",0.25);
                dataVec[i][2] = parameterList->sublist("Parameter").get("initial volume fraction solid material"+std::to_string(i+1),0.5); //do we need this?
                dataVec[i][3] = parameterList->sublist("Parameter").get("Darcy parameter material"+std::to_string(i+1),1.e-2);
                dataVec[i][4] = parameterList->sublist("Parameter").get("Youngs modulus material"+std::to_string(i+1),60.e6);
                dataVec[i][5] = parameterList->sublist("Parameter").get("Poisson ratio material"+std::to_string(i+1),0.3);
            }
            else if (dim==3) {
                dataVec[i].resize(12);
                dataVec[i][0] = parameterList->sublist("Parameter").get("Youngs modulus material"+std::to_string(i+1),2.e5);
                dataVec[i][1] = parameterList->sublist("Parameter").get("Poisson ratio material"+std::to_string(i+1),0.3);
                dataVec[i][2] = 0.; //body force x
                dataVec[i][3] = 0.; //body force y
                dataVec[i][4] = parameterList->sublist("Parameter").get("body force z"+std::to_string(i+1),0.);; //body force z
                dataVec[i][5] = parameterList->sublist("Parameter").get("initial volume fraction solid material"+std::to_string(i+1),0.67);
                dataVec[i][6] = parameterList->sublist("Parameter").get("Darcy parameter material"+std::to_string(i+1),0.01);
                dataVec[i][7] = 2000.; //effective density solid
                dataVec[i][8] = 1000.; //effective density fluid?
                dataVec[i][9] = 9.81;  // gravity
                dataVec[i][10] = parameterList->sublist("Timestepping Parameter").get("Newmark gamma",0.5);
                dataVec[i][11] = parameterList->sublist("Timestepping Parameter").get("Newmark beta",0.25);
            }
        }
        
        
        else if (tpmType == "Biot-StVK") {
            dataVec[i][0] = parameterList->sublist("Parameter").get("Youngs modulus material"+std::to_string(i+1),60.e6);
            dataVec[i][1] = parameterList->sublist("Parameter").get("Poisson ratio material"+std::to_string(i+1),0.3);
            dataVec[i][2] = parameterList->sublist("Parameter").get("initial volume fraction solid material"+std::to_string(i+1),0.5); //do we need this?
            dataVec[i][3] = parameterList->sublist("Parameter").get("Darcy parameter material"+std::to_string(i+1),1.e-2);
            dataVec[i][4] = parameterList->sublist("Timestepping Parameter").get("Newmark gamma",0.5);
            dataVec[i][5] = parameterList->sublist("Timestepping Parameter").get("Newmark beta",0.25);
        }
    }
    ///---------------------------------------------------------------------------------------------------------------------------------------
    ///---------------------------------------------------------------------------------------------------------------------------------------
    //              fuellen der element specs
    ///---------------------------------------------------------------------------------------------------------------------------------------
    for (int i=0; i<es_vec.size(); i++){
        if(tpmType == "Biot"){
            if (dim==2)
                this->SMTSetElSpecBiot(&es_vec[i] ,&idata, ic, ng, dataVec[i]);
            else if(dim==3)
                this->SMTSetElSpecBiot3D(&es_vec[i] ,&idata, ic, ng, dataVec[i]);
        }
        else if(tpmType == "Biot-StVK")
            this->SMTSetElSpecBiotStVK(&es_vec[i] ,&idata, ic, ng, dataVec[i]);
    }

    ///---------------------------------------------------------------------------------------------------------------------------------------
    //              processing variables - eigentlich sollte das hier aus elmt spec sein
    ///---------------------------------------------------------------------------------------------------------------------------------------
    LO elementSizePhase = elements1->nodesPerElement();
    LO sizePhase = dim * elementSizePhase;
    LO sizePressure = elements2->nodesPerElement();
    GO sizePhaseGlobal = A00->getMap()->getMaxAllGlobalIndex()+1;

    ///---------------------------------------------------------------------------------------------------------------------------------------
    //              prepare working vector
    ///---------------------------------------------------------------------------------------------------------------------------------------
    int workingVectorSize;
    if(tpmType == "Biot"){
        if (dim==2)
            workingVectorSize = 5523;
        else if(dim==3)
            workingVectorSize = 1817;
    }
    else if(tpmType == "Biot-StVK")
        workingVectorSize = 5223;
    
    double* v = new double [workingVectorSize];

    // nd sind Nodalwerte, Anzahl an structs in nd sollte den Knoten entsprechen, bei P2-P1 in 2D also 9
    // In X stehen die Koordinaten, X[0] ist x-Koordinate, X[1] ist y-Koordinate, etc.
    // nd->X[0]
    // at ist die Loesung im letzten Newtonschritt.
    // nd[0]->at[0];
    // ap ist die Loesung im letzten Zeitschritt.
    // nd[0]->ap[0]
    // rdata ist die Zeitschrittweite, RD_TimeIncrement wird in sms.h definiert, entsprechend wird auch die Laenge von rdata dort definiert. Standard 400, aber auch nicht gesetzt. Wert muss selber initialisiert werden; eventuell kuerzer moeglich.
    ///---------------------------------------------------------------------------------------------------------------------------------------
    //              erstellen von rdata  ! sieht aus als waere in der feddlib nur constante zeitschrittweite erlaubt
    ///---------------------------------------------------------------------------------------------------------------------------------------
    std::vector<double> rdata(RD_TimeIncrement+1, 0.);

    rdata[RD_TimeIncrement] = parameterList->sublist("Timestepping Parameter").get("dt",0.01);

    ///---------------------------------------------------------------------------------------------------------------------------------------
    //              restelen von node spec - hier nur dummy pointer 
    ///---------------------------------------------------------------------------------------------------------------------------------------
    NodeSpec *ns=NULL;//dummy not need in SKR

    ///---------------------------------------------------------------------------------------------------------------------------------------
    //              erstellen eines node data vectors ! warum hier eigentlich dann nciht als vector ...
    ///---------------------------------------------------------------------------------------------------------------------------------------
        
    NodeData** nd = new NodeData*[ elementSizePhase + sizePressure ];

    for (int i=0; i<elementSizePhase + sizePressure; i++){
        nd[i] = new NodeData();
    }

    int numNodes = elementSizePhase + sizePressure;
    
    // vec2D_dbl_Type ==> std::vector<std::vector<double> >
    // hier werden matrizen erstellen e.g. XI(NoNodes, Dim) und mit null initialisiert
    // diese stellen den permanenten speicher da, danach wird die erstelte node data
    // auf diesen speicher verwiesen ! kann man hier nicht direkt in die feddlib verweisen
    vec2D_dbl_Type xFull( numNodes, vec_dbl_Type(dim,0.) );
    vec2D_dbl_Type atFull( numNodes, vec_dbl_Type(dim,0.) );
    vec2D_dbl_Type apFull( numNodes, vec_dbl_Type(dim,0.) );
    
    for (int i=0; i<elementSizePhase + sizePressure; i++) {
        nd[i]->X = &(xFull[i][0]);
        nd[i]->at = &(atFull[i][0]);
        nd[i]->ap = &(apFull[i][0]);
    }

    //das hier sind die kooridinaten ... sortiert wie freiheitsgrade etwas weiter unten
    GO offsetMap1 = dim * mapRepeated1->getMaxAllGlobalIndex()+1;
    vec2D_dbl_ptr_Type pointsRepU = domainVec_.at(0)->getPointsRepeated();
    vec2D_dbl_ptr_Type pointsRepP = domainVec_.at(1)->getPointsRepeated();
    
    // das hier isnd freiheitsgrade uArrayNewton aktuell (at)
    // und alte freiheitsgrade      uArrayTime           (ap)
    // speicherreihenfolge ist: uArrayTime = [ n1x, n1y, n1z, ... nNx, nNy, nNz]
    // for N -> global node number !!
    Teuchos::ArrayRCP< const SC > uArrayNewton = u_repeatedNewton->getData(0);
    Teuchos::ArrayRCP< const SC > pArrayNewton = p_repeatedNewton->getData(0);
    Teuchos::ArrayRCP< const SC > uArrayTime = u_repeatedTime->getData(0);
    Teuchos::ArrayRCP< const SC > pArrayTime = p_repeatedTime->getData(0);

    // aufbauen der localen steifigkeitsmatrix 
    double** mat = new double*[sizePhase+sizePressure];
    for (int i=0; i<sizePhase+sizePressure; i++){
        mat[i] = new double[sizePhase+sizePressure];
    }
    
    // das sind die rechten seiten!
    // diese werden allerdings nicht mehr veraendert sondern der pointer unten 
    // wird zur modifikation verwendet!!
    Teuchos::ArrayRCP<SC> fValues0 = F0->getDataNonConst(0);
    Teuchos::ArrayRCP<SC> fValues1 = F1->getDataNonConst(0);
    
    // Element loop

    ElementData ed = ElementData();
    for (UN T=0; T<elements1->numberElements(); T++) {
        
        std::vector<double> tmpHp = ed_[T]->getHp(); // Dies sind die alten Daten
        std::vector<double> tmpHt = ed_[T]->getHt(); // Dies sind die neuen Daten nachdem das Element aufgerufen wurde, wir hier eigentlich nicht als Variable in ed_ benoetigt.
        ed.hp = &tmpHp[0];
        ed.ht = &tmpHt[0];
        
        int materialFlag = elements1->getElement(T).getFlag();
        TEUCHOS_TEST_FOR_EXCEPTION( materialFlag>es_vec.size()-1, std::runtime_error, "There are not enought material parameters initialized." ) ;
        int counter=0;
        //Newtonloesung at und Zeitschrittloesung ap
        for (int j=0; j<elementSizePhase; j++) {
            for (int d=0; d<dim; d++) {
                ////////-----------------------------
                // Interesting, hiere is the mesh!!! ==> elements1->getElement(T).getNode(j)
                ////////-----------------------------
                LO index = dim * elements1->getElement(T).getNode(j)+d;//dim * elements1->at(T).at( j ) + d;
                atFull[j][d] = uArrayNewton[index];
                apFull[j][d] = uArrayTime[index];
            }
        }
        for (int j=0; j<sizePressure; j++) {
            LO index = elements2->getElement(T).getNode(j);//elements2->at(T).at( j );
            atFull[elementSizePhase+j][0] = pArrayNewton[index];
            apFull[elementSizePhase+j][0] = pArrayTime[index];
        }
        
        //Nodes
        for (int j=0; j<elementSizePhase; j++ ) {
            LO index = elements1->getElement(T).getNode(j);
            for (int d=0; d<dim; d++) {
                xFull[j][d] = (*pointsRepU)[index][d];
            }
        }
        for (int j=0; j<sizePressure; j++ ) {
            LO index = elements2->getElement(T).getNode(j);
            for (int d=0; d<dim; d++) {
                xFull[elementSizePhase+j][d] = (*pointsRepP)[index][d];
            }
        }
        // ab hier ist alles in node data definiert (allerdings nur at, ap und X)

        // allocieren des localen element vectors
        vec_dbl_Type p( sizePhase+sizePressure , 0. );

        // initialization of local element matrix 
        for (int i=0; i<sizePhase+sizePressure; i++){
            for (int j=0; j<sizePhase+sizePressure; j++)
                mat[i][j] = 0.;
        }

        // hier kann nun auch ein anderes element stehen !!!
        // element assembly
        if(tpmType == "Biot"){
            if(dim==2)
                this->SKR_Biot( v, &es_vec[materialFlag], &ed, &ns, nd , &rdata[0], &idata, &p[0], mat  );
            else if (dim==3)
                this->SKR_Biot3D( v, &es_vec[materialFlag], &ed, &ns, nd , &rdata[0], &idata, &p[0], mat  );
        }
        else if(tpmType == "Biot-StVK")
            this->SKR_Biot_StVK( v, &es_vec[materialFlag], &ed, &ns, nd , &rdata[0], &idata, &p[0], mat  );
        

        // history update need to be extended for rdata, idata and nodal history
        if (updateHistory)
            ed_[T]->setHp( ed.ht );
        

    ///---------------------------------------------------------------------------------------------------------------------------------------
    //              ab hier werden locale matrix und vector assembliert (im localen rank!)
    ///---------------------------------------------------------------------------------------------------------------------------------------
    
    // sizePhase  -> anzahl element dofs aus dem ersten feld (dsip)
    // sizePressure -> anzahl der element dofs fuer zweites feld (druck)
    // globale steifigkeit vermutlich so gedacht fuer zwei felder
    //  
    //              A = |A00 A01|
    //                  |A10 A11|
    //
        if (update) {
                    
            // A00 & A01
            for (UN i=0; i < sizePhase; i++) {
                Teuchos::Array<SC> value00( sizePhase, 0. );
                Teuchos::Array<GO> indices00( sizePhase, 0 );
                for (UN j=0; j < value00.size(); j++) {
                    
                    value00[j] = mat[i][j];
                    
                    LO tmpJ = j/dim;
                    LO index = elements1->getElement(T).getNode(tmpJ);
                    if (j%dim==0)
                        indices00[j] = dim * mapRepeated1->getGlobalElement( index );
                    else if (j%dim==1)
                        indices00[j] = dim * mapRepeated1->getGlobalElement( index ) + 1;
                    else if (j%dim==2)
                        indices00[j] = dim * mapRepeated1->getGlobalElement( index ) + 2;
                }
                
                Teuchos::Array<SC> value01( sizePressure, 0. );
                Teuchos::Array<GO> indices01( sizePressure, 0 );

                for (UN j=0; j < value01.size(); j++) {
                    value01[j] = mat[i][sizePhase+j];
                    LO index = elements2->getElement(T).getNode(j);
                    indices01[j] = mapRepeated2->getGlobalElement( index );
                }
                
                GO row;
                LO tmpI = i/dim;
                LO index = elements1->getElement(T).getNode(tmpI);
                if (i%dim==0)
                    row = dim * mapRepeated1->getGlobalElement( index );
                else if (i%dim==1)
                    row = dim * mapRepeated1->getGlobalElement( index ) + 1;
                else if (i%dim==2)
                    row = dim * mapRepeated1->getGlobalElement( index ) + 2;
                
                A00->insertGlobalValues( row, indices00(), value00() );
                A01->insertGlobalValues( row, indices01(), value01() );
                
                if (i%dim==0)
                    fValues0[ dim*index ] += p[ i ];
                else if (i%dim==1)
                    fValues0[ dim*index+1 ] += p[ i ];
                else if (i%dim==2)
                    fValues0[ dim*index+2 ] += p[ i ];
            }
            // A10 & A11
            for (UN i=0; i < sizePressure; i++) {
                Teuchos::Array<SC> value10( sizePhase   , 0. );
                Teuchos::Array<GO> indices10( sizePhase   , 0 );
                for (UN j=0; j < value10.size(); j++) {
                    value10[j] = mat[sizePhase+i][j];
                    
                    LO tmpJ = j/dim;
                    LO index = elements1->getElement(T).getNode(tmpJ);
                    if (j%dim==0)
                        indices10[j] = dim * mapRepeated1->getGlobalElement( index );
                    else if (j%dim==1)
                        indices10[j] = dim * mapRepeated1->getGlobalElement( index ) + 1;
                    else if (j%dim==2)
                        indices10[j] = dim * mapRepeated1->getGlobalElement( index ) + 2;
                }
                
                Teuchos::Array<SC> value11( sizePressure, 0. );
                Teuchos::Array<GO> indices11( sizePressure, 0 );
                for (UN j=0; j < value11.size(); j++) {
                    value11[j] = mat[sizePhase+i][sizePhase+j];
                    
                    LO index = elements2->getElement(T).getNode(j);
                    indices11[j] = mapRepeated2->getGlobalElement( index );
                }

                
                LO index2 = elements2->getElement(T).getNode(i);
                GO row = mapRepeated2->getGlobalElement( index2 );
                A10->insertGlobalValues( row, indices10(), value10() );
                A11->insertGlobalValues( row, indices11(), value11() );
                
                fValues1[ index2 ] += p[ sizePhase + i ];
            }
        }
    }
    
    // free intermediate vairable memory
    for (int i=0; i<sizePhase+sizePressure; i++)
        delete [] mat[i];
    delete [] mat;
    
    delete [] v;
    
    for (int i=0; i<elementSizePhase+sizePressure; i++)
        delete nd[i];
    
    delete [] nd;
    
    // das sind bestimmt die parallel synchonisations befehle
    // heist das eigentlich dass man das fuer einen vector nicht machen muss ?
    A00->fillComplete( A00->getMap("row"), A00->getMap("row") );
    A01->fillComplete( A10->getMap("row"), A00->getMap("row") );
    A10->fillComplete( A00->getMap("row"), A10->getMap("row") );
    A11->fillComplete( A10->getMap("row"), A10->getMap("row") );
    
}