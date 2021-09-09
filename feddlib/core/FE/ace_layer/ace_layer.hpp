#pragma once
// Hard Code path otherwise configured by cmake
#define SMSHeader //anpassen!!
#define SMSUtil //anpassen!!
#define SMSCCompiler "gcc"

/*
A wrapper for using AceGen elements in FEDDlib.
Notice that in this context an object of type AceGenElement
does not represent a single element in a discretization but 
just holds information and gives access to the selected elements
functions. A similar concept is the Domain concept in 
AceFEM defined via SMTAddDomain.
The concept for an AceGen/AceFEM user is to build elements
with the "Environment" -> "AceFEM" setting, with as 
few limitations as possible.
To ensure inviolability of the ElementSpecification (which 
are FIXED by the element code) the AceGenElement must use
getter functions for securely exchange information with the rest of the code.
*/

#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <stdlib.h>

extern "C"
{
    #include "sms.h"
}

// in view of MPI parallelism, std::cout is only used when building in debug config (cmake)
#ifndef CMakeDebugConfig
    #define CMakeDebugConfig 0
#endif

typedef int (CALLBACK* DLLSetElSpec)(ElementSpec *,int *,int i,int j);

// Interface Class for convenient access of an AceFEM element
// provided by a AceGen generated c-code.
class AceGenElement
{
    public:
        AceGenElement() : es(), shared_elmt_lib(nullptr), compiler_call("gcc -shared -fpic") {};
        ~AceGenElement();
        AceGenElement(const AceGenElement&) = delete;
        AceGenElement& operator=(const AceGenElement&) = delete;

    public:
        // Constructor with automatic initialization.
        // Input is the path to the element c-code.
        // Example: AceGenElement("~/my_element.c");
        AceGenElement(std::string path_to_elememt);
    
    public:
        // Load an AceGen generated element.
        // The input MUST be the parth to the c-code.
        // At load(), the code will be compiled always.
        // After succesfull load(), this AceGenElement has
        // knowledge about the AceGen elements properties,
        // e.g., number of nodes or length of histroy,
        // and access to its defined functions, e.g., SKR().
        int load(std::string path_to_elememt);

    public:
        // Releases a loaded shared element library.
        int unload();

    public:
        // Prints a summary of all knwon properties of this AceGenElement.
        int status();

    public:
        // AceGenElement initialization function.
        // The element is loaded from a shared library and
        // does other initialization.
        int init(std::string path_to_elememt);

    public:
        // Calls the SKR - Tangent and Residual function of this element type.
        // Writes the element stiffness matrix and load vector to the s and p input.
        // The load vector is ordered InternalForceVector = [p1, p2, ..., px] 
        // for x elemental d.o.f..
        // The stiffness matrix is exported in a flattened form:
        // StiffnessVector = [p11, p12, ..., p1x, p21, p22, ..., p2x, ..., pyx] 
        // for x and y elemental d.o.f..
        // Returns 0;
        // Hint: Access is s[2][3] = StiffnessVector[2*y+3]
        // Hint: Initialization of Stiffness- and InternalForceVector not required.
        int SKR(ElementData* Element_Data, NodeData** VectorOfNodeData,
                double* SystemRealData, int* SystemIntegerData,
                double* InternalForceVector, double* StiffnessVector
                );

    private:
        // Internal Data structure for element specifications.
        ElementSpec es;

    private:
        // Pointer to the currently loaded shared element library.
        void* shared_elmt_lib;

    public:
        // String to be executed in a shell to compile a shared library.
        // Default: "gcc -shared -fpic"
        std::string compiler_call;

    private:
        // AceGen working vector, allocated once for this type of element.
        double* v;

    public:
        // ElementSpecification Interface:
        int get_NoNodes() const {return es.id.NoNodes;}
        int get_NoDofs() const {return es.id.NoDOFGlobal;}
        int get_NoNodalDofs(int i) const {return this->es.DOFGlobal[i];};
        double* get_ReferenceNodes() const {return &(es.ReferenceNodes[0]);}
        size_t get_NoTimeStorage() const {return this->es.id.NoTimeStorage;}
        size_t get_NoElementData() const {return this->es.id.NoElementData;}
};

// Represents the AceFEM datastructures such as nodedata
// elementdata idata and rdata
// provides update routines according to expected situations
// e.g. update_step -> switches history etc
//      update_iteration -> updates current nodal solution etc
// in a parallel setup this can be hold per rank
/* 
Notation Remark:
        ap -> previous value of degrees of freedom (equilibrium at the end of last step)
        at -> current value of degrees of freedom (at = ap + da, at the beginning of the time step at=ap)

Data that need to be updated per time step/ step:
    Nodal dof fields                ap <- at
    Elment/Nodal history fields     hp <- ht
    Current BC multiplier           Bp <- Bt    (later)
    Nodal Sensitivity Data          sp <- st    (later)


Data that need to be updated after each linear solution (e.g. iteration from Newton scheme)
    Nodal dof fields                da <- da
    Nodal dof fields                at <- ap + da

Important general data:
    idata$$["NoDimensions"]         number of spatial dimensions of the problem (2 or 3)
    idata$$["Iteration"]            index of the current iteration within the iterative loop
    idata$$["Step"]                 total number of completed solution steps (set by Newton-Raphson iterative procedure)
    idata$$["ErrorStatus"]          look for errors here
    idata$$["SubDivergence"]        look for errors here
    
    rdata$$["Multiplier"]           current global load multiplier (relevant within the elements e.g. self weight)
    rdata$$["Time"]                 current real time for which an equilibrium is solved for
    rdata$$["TimeIncrement"]        current time increment from the last knwon equilibrium to the real time for which an equilibrium is solved for

Corresponding codes:
    
    #define ID_NoDimensions         24
    #define ID_Iteration            5
    #define ID_Step                 63
    #define ID_ErrorStatus          8
    #define ID_SubDivergence        15

    #define RD_Multiplier            0
    #define RD_Time                  5
    #define RD_TimeIncrement         6
*/


// This class handels the data for elements.
// Those are fields hp,ht and Data;
class AceFEMElementDataContainer
{
    public:
        AceFEMElementDataContainer(size_t NoElements, ElementSpec* ElementSpecifications) :
        ElementDataContainer(NoElements), 
        ht_Container(NoElements*(ElementSpecifications->id.NoTimeStorage), 0.0),
        hp_Container(NoElements*(ElementSpecifications->id.NoTimeStorage), 0.0),
        Data_Container(NoElements*(ElementSpecifications->id.NoElementData), 0.0)
        {
            size_t NoTS(ElementSpecifications->id.NoTimeStorage);
            size_t NoED(ElementSpecifications->id.NoElementData);
            for (size_t i=0; i<NoElements; i++)
            {
                ElementDataContainer[i].ht = &ht_Container[i*NoTS];
                ElementDataContainer[i].hp = &hp_Container[i*NoTS];
                ElementDataContainer[i].Data = &Data_Container[i*NoED];
            }
        };

    public:
        // Update time dependent values: hp <- ht .
        int time_update()
        {
            hp_Container.swap(ht_Container);
            return 0;
        };

    private:
        std::vector<ElementData> ElementDataContainer;
        std::vector<double> ht_Container;
        std::vector<double> hp_Container;
        std::vector<double> Data_Container;
};

// This class handels the data for nodes.
// Those are fields hp,ht and Data;
class AceFEMNodeDataContainer
{
    public:
        AceFEMNodeDataContainer(size_t NoNodes, NodeSpec* NodeSpecifications) :
        NodeDataContainer(NoNodes),
        hp_Container(NoNodes*(NodeSpecifications->id.NoNodeStorage), 0.0),
        ht_Container(NoNodes*(NodeSpecifications->id.NoNodeStorage), 0.0),
        Data_Container(NoNodes*(NodeSpecifications->id.NoData), 0.0)
        {
            size_t NoTS(NodeSpecifications->id.NoNodeStorage);
            size_t NoND(NodeSpecifications->id.NoNodeData);
            for (size_t i=0; i<NoNodes; i++)
            {
                NodeDataContainer[i].ht = &ht_Container[i*NoTS];
                NodeDataContainer[i].hp = &hp_Container[i*NoTS];
                NodeDataContainer[i].Data = &Data_Container[i*NoND];
            }
        };

    public:
        // Update time dependent values: hp <- ht .
        int time_update()
        {
            hp_Container.swap(ht_Container);
            return 0;
        };

    private:
        std::vector<NodeData> NodeDataContainer;
        std::vector<double> ht_Container;
        std::vector<double> hp_Container;
        std::vector<double> Data_Container;
};

// Management layer of Data that AceFEM elements expect.
class AceFEMDataLayer
{
    public:
        AceFEMDataLayer(size_t NoElements, size_t NoNodes){};

    public:
        //SystemRealData
        double SystemRealData[10];
        //SystemRealData
        int SystemIntegerData[100];

    public:
        // Updates all time dependent data of Nodes and Elements.
        // Also updates required Real and Integer data.
        int time_step(double CurrentTimeIncrement, double CurrentLoadMultiplier);


};