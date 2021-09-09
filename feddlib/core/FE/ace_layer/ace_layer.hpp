#pragma once
// Hard Code path otherwise configured by cmake
#define SMSHeader "/Users/cnisters/source/feddlib_mech/feddlib/core/FE/ace_layer/"
#define SMSUtil "/Users/cnisters/source/feddlib_mech/feddlib/core/FE/ace_layer/SMSUtility.c"
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
#include <iomanip>
#include <string>
#include <vector>
#include <array>
#include <stdlib.h>

// #include "sms.h"

// extern "C"
// {
//     #include "sms.h"
// }

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

inline int AceFEMDataLayer::time_step(double CurrentTimeIncrement, double CurrentLoadMultiplier)
{
    // update step dependent integer data
    this->SystemIntegerData[ID_Step] += 1;
    this->SystemIntegerData[ID_Iteration] = 0;

    // update step dependent real data
    this->SystemRealData[RD_TimeIncrement] = CurrentTimeIncrement;
    this->SystemRealData[RD_Time] += CurrentTimeIncrement;
    this->SystemRealData[RD_Multiplier] = CurrentLoadMultiplier;

    return 0;
};


inline int AceGenElement::load(std::string path_to_elememt)
{
    // if a shared library is loaded, it must be released frist
    if (shared_elmt_lib!=nullptr)
    {
        std::cout << "Warning: Unloading previously loaded element!" << std::endl;
        int unld = this -> unload();
        if (unld!=0)
        {
            std::cerr << "Error: Failed to unload a previously loaded element!" << std::endl;
            return -5;
        }
    }

    if (CMakeDebugConfig) std::cout << "Loading Element: " << path_to_elememt << std::endl;
    
    // always recompile the element
    std::size_t dir = path_to_elememt.rfind("/");
    if (dir==std::string::npos)
    {
        std::cerr << "Error: Element path: " << path_to_elememt << " is not a valid path." << std::endl;
        return -1;
    }
    //compiler_call = SMSCCompiler; // future feature set by cmake
    // and stop copying ! instead we give the compiler the paths to the assets !
    std::string element_directory(path_to_elememt.substr(0, dir+1));
    std::string element_name(path_to_elememt.substr(dir+1, path_to_elememt.size()));
    if (CMakeDebugConfig) std::cout << "Element Directory: " << element_directory << std::endl;
    if (CMakeDebugConfig) std::cout << "Element Name     : " << element_name << std::endl;
    std::string path_to_elememt_shared_object(path_to_elememt);
    path_to_elememt_shared_object.pop_back();
    path_to_elememt_shared_object.append("so");
    if (CMakeDebugConfig) std::cout << "Compiling Element: " << path_to_elememt_shared_object << std::endl;
    std::string compile_command(compiler_call + " -o" + path_to_elememt_shared_object + " " + path_to_elememt+ " " + SMSUtil + " -I" + SMSHeader);
    if (CMakeDebugConfig) std::cout << "Command: " << compile_command << std::endl;
    int ret = system(&compile_command[0]);
    if (ret==-1)
    {
        std::cerr << "Error: Compilation of " << path_to_elememt << " failed!" << std::endl;
        return -2;
    }

    // loading the newly generated shared library
    char* error;
    shared_elmt_lib = dlopen(&path_to_elememt_shared_object[0], RTLD_LAZY);
    error = dlerror ();
    if (error)
    {
        std::cerr << "Error: Loading shared element library " << path_to_elememt_shared_object << " failed!" << std::endl;
        return -3;    
    }
    DLLSetElSpec loaded_element_setup=(DLLSetElSpec)dlsym(shared_elmt_lib,"SMTSetElSpec");
    error = dlerror ();
    if (error)
    {
        std::cerr << "Error: Loading functions from shared element library " << path_to_elememt_shared_object << " failed!" << std::endl;
        return -4;    
    }
    this->es.user.SetElSpec=loaded_element_setup;

    // int SMTSetElSpec(ElementSpec *es,int *idata,int ic,int ng)
    // This function sets up the ElementSpec and returns 0.
    // *idata -> idata[8] is just filled with an error code (3) if unvalid integration code. nullptr argument is possible
    // ic -> possible integration code for overwriting the one defined in the element. ic=-1 uses the element integration code.
    // ng -> possible default domain data for overwriting the one defined in the element. ng=-1 uses the element domain data.
    // Remarks: Domain data modification (ng other than -1) leaves the es->Data=nullptr. So do this ONLY if you know what you do.
    //          This function sets up properly, but integration points must be assigned in the following!
    
    int idata[9]; // Dummy idata here ... might be changed later!
    this->es.user.SetElSpec(&(this->es), idata, -1, -1);
    this->es.IntPoints = SMTMultiIntPoints(&(this->es.id.IntCode), idata, &(this->es.id.NoIntPoints), &(this->es.id.NoIntPointsA), &(this->es.id.NoIntPointsB), &(this->es.id.NoIntPointsC), 1);
    if (idata[ID_ErrorStatus]==ERROR_Fatal)
    {
        std::cerr << "Error: Processing AceFEM integration points failed! Possible issue with unsupported integration codes. See: [AceGen/tutorial/Numerical Integration]" << std::endl;
        return -6;    
    }

    // Some infos about the integration functions from SMSUtility.c

    // double *SMTIntPoints(int *icode,int *ngpo)
    // Returns a pointer to the vector of integration points.
    // icode -> actual input rendering the integration code
    // ngpo -> is modified towads number of integration points

    // double *SMTMultiIntPoints(int *icode,int *idata,int *ngpo, int *na,int *nb,int *nc,int alloc)
    // Returns a pointer to the vector of integration points (if alloc = 1!!!).
    // alloc -> int for wether a new double vector (which is filled by this function with gauss points) has to be allocated
    // *idata -> idata[8] is just filled with an error code (3) if unvalid integration code. nullptr argument is possible
    // ngpo -> is modified towads resulting number of integration points
    // icode -> actual input rendering the integration code - if NoIntPointsA - NoIntPointsB - NoIntPointsC are used the integartion code is in [202020, ..., 343434]
    // na, nb, nc -> will be modified to NoIntPointsA - NoIntPointsB - NoIntPointsC

    return 0;
};

inline int AceGenElement::SKR(
        ElementData* Element_Data, NodeData** VectorOfNodeData,
        double* SystemRealData, int* SystemIntegerData,
        double* InternalForceVector, double* StiffnessVector
        )
{
    int no_dofs(this->es.id.NoDOFGlobal);
    int no_nodes(this->es.id.NoNodes);
    
    // initialize output fields and temporary matrix pointers
    double * s_[no_dofs];
    for (int i=0; i<no_dofs; i++)
    {
        InternalForceVector[i] = 0.0;
        for (int j=0; j<no_dofs; j++) StiffnessVector[i*no_dofs+j] = 0.0;
        s_[i] = &StiffnessVector[i*no_dofs];
    }

    // Some notes on the AceFEM SKR routine:
    // v -> workingVector (allocated by this AceGenElement)
    // es -> element specification of this element type (this AceGenElement specifically)
    // ed -> element data for the current element in a discretization for which to evaluate SKR
    // ns -> node specification (only barely used in AceFEM and its functionality is covered by FEDDlib for the most part HENCE NOT SUPPORTED BY NOW)
    // nd -> vector of pointer to NodeData for each node of this element.
    // rdata -> vector of global simulation double parameters, e.g., current time
    // idata -> vector of global simulation integer parameters, e.g., current step
    // p -> local element vector / internal force vector
    // s -> loacal element matrix / stiffness matrix
    // SKR(double v[110],ElementSpec *es,ElementData *ed,NodeSpec **ns,NodeData **nd,double *rdata,int *idata,double *p,double **s)


    // call the element SKR routine
    this->es.user.SKR(this->v, &(this->es), Element_Data, nullptr, VectorOfNodeData, SystemRealData, SystemIntegerData, &InternalForceVector[0], &s_[0]);

    // handle symmetic matrices (as long as not taken advantage of from FEDDlib)
    if (this->es.id.SymmetricTangent)
        for (int i=0; i<no_dofs; i++)
            for (int j=i; j<no_dofs; j++) 
                StiffnessVector[i*no_dofs+j] = StiffnessVector[j*no_dofs+i];


    //debug putput
    for (int i=0; i<no_dofs; i++)
    {
        for (int j=0; j<no_dofs; j++) 
                if (CMakeDebugConfig) std::cout << " " << StiffnessVector[i*no_dofs+j];
        if (CMakeDebugConfig) std::cout << std::endl;
    }
            

    return 0;
};

inline int AceGenElement::unload()
{
    if (shared_elmt_lib!=nullptr) dlclose(shared_elmt_lib);  
    if (CMakeDebugConfig) std::cout << "AceGenElement unloaded!" << std::endl;
    return 0;
};

inline AceGenElement::AceGenElement(std::string path_to_elememt) : AceGenElement()
{
    init(path_to_elememt);
};

inline AceGenElement::~AceGenElement()
{
    unload();
    delete this->v;
};

inline int AceGenElement::init(std::string path_to_elememt)
{
    load(path_to_elememt);
    this -> v = new double [this->es.id.WorkingVectorSize];
    return 0;
};

inline int AceGenElement::status()
{
    std::cout << "Element Status: " << std::endl;
    std::cout << "\t SpecIndex:               " << this->es.id.SpecIndex << std::endl;
    std::cout << "\t NoDimensions:            " << this->es.id.NoDimensions << std::endl;
    std::cout << "\t NoDOFGlobal:             " << this->es.id.NoDOFGlobal << std::endl;
    std::cout << "\t NoDOFCondense:           " << this->es.id.NoDOFCondense << std::endl;
    std::cout << "\t NoNodes:                 " << this->es.id.NoNodes << std::endl;
    std::cout << "\t DOFGlobal:               {";
    for (int i=0; i<(this->es.id.NoNodes); i++) std::cout << " " << this->es.DOFGlobal[i];
    std::cout << " }" << std::endl;
    std::cout << "\t NoDomainData:            " << this->es.id.NoDomainData << std::endl;
    std::cout << "\t NoSegmentPoints:         " << this->es.id.NoSegmentPoints << std::endl;
    std::cout << "\t IntCode:                 " << this->es.id.IntCode << std::endl;
    std::cout << "\t NoTimeStorage:           " << this->es.id.NoTimeStorage << std::endl;
    std::cout << "\t NoElementData:           " << this->es.id.NoElementData << std::endl;
    std::cout << "\t NoGPostData:             " << this->es.id.NoGPostData << std::endl;
    std::cout << "\t NoNPostData:             " << this->es.id.NoNPostData << std::endl;
    std::cout << "\t SymmetricTangent:        " << this->es.id.SymmetricTangent << std::endl;
    std::cout << "\t NoSensNames:             " << this->es.id.NoSensNames << std::endl;
    std::cout << "\t ShapeSensitivity:        " << this->es.id.ShapeSensitivity << std::endl;
    std::cout << "\t NoIData:                 " << this->es.id.NoIData << std::endl;
    std::cout << "\t NoRData:                 " << this->es.id.NoRData << std::endl;
    std::cout << "\t DefaultIntegrationCode:  " << this->es.id.DefaultIntegrationCode << std::endl;
    //std::cout << "\t LocalReKe:               " << this->es.id.LocalReKe << std::endl;
    std::cout << "\t NoAdditionalData:        " << this->es.id.NoAdditionalData << std::endl;
    std::cout << "\t NoCharSwitch:            " << this->es.id.NoCharSwitch << std::endl;
    std::cout << "\t NoIntSwitch:             " << this->es.id.NoIntSwitch << std::endl;
    std::cout << "\t NoDoubleSwitch:          " << this->es.id.NoDoubleSwitch << std::endl;
    std::cout << "\t PostIterationCall:       " << this->es.id.PostIterationCall << std::endl;
    std::cout << "\t Active:                  " << this->es.id.Active << std::endl;
    std::cout << "\t DOFScaling:              " << this->es.id.DOFScaling << std::endl;
    std::cout << "\t EBCSensitivity:          " << this->es.id.EBCSensitivity << std::endl;
    std::cout << "\t SensitivityOrder:        " << this->es.id.SensitivityOrder << std::endl;
    std::cout << "\t WorkingVectorSize:       " << this->es.id.WorkingVectorSize << std::endl;
    std::cout << "\t Topology:                " << std::string(this->es.Topology) << std::endl;
    std::cout << "\t DomainDataNames:         " << std::endl;
    for (int i=0; i<(this->es.id.NoDomainData); i++)
        std::cout << "\t\t\t\t\t\t\t  - " << "\"" << std::string(this->es.DomainDataNames[i]) << "\" -> " << this->es.Data[i] << std::endl;
    std::cout << "\t GPostNames:              " << std::endl;
    for (int i=0; i<(this->es.id.NoGPostData); i++)
        std::cout << "\t\t\t\t\t\t\t  - " << "\"" << std::string(this->es.GPostNames[i]) << "\"" << std::endl;
    std::cout << "\t NPostNames:              " << std::endl;
    for (int i=0; i<(this->es.id.NoNPostData); i++)
        std::cout << "\t\t\t\t\t\t\t  - " <<"\"" << std::string(this->es.NPostNames[i]) << "\"" << std::endl;
    std::cout << "\t Segments:                " << *(this->es.Segments) << std::endl;
    std::cout << "\t NoNodeStorage:           " << *(this->es.NoNodeStorage) << std::endl;
    std::cout << "\t NoNodeData:              " << *(this->es.NoNodeData) << std::endl;
    std::cout << "\t NodeID:                  " << std::endl;
    for (int i=0; i<(this->es.id.NoNodes); i++)
        std::cout << "\t\t\t\t\t\t\t  - " << "\"" << std::string(this->es.NodeID[i]) << "\"" << std::endl;
    std::cout << "\t ReferenceNodes:          " << std::endl;
    for (int i=0; i<(this->es.id.NoNodes); i++)
        std::cout << "\t\t\t\t\t\t\t  X_"<< i << ":\t [ " << this->es.ReferenceNodes[0+i*3] << ",\t" << this->es.ReferenceNodes[1+i*3] << ",\t" << this->es.ReferenceNodes[2+i*3] << " ]" << std::endl;
    std::cout << "\t NoIntPoints:             " << this->es.id.NoIntPoints << std::endl;
    std::cout << "\t NoIntPointsA:            " << this->es.id.NoIntPointsA << std::endl;
    std::cout << "\t NoIntPointsB:            " << this->es.id.NoIntPointsB << std::endl;
    std::cout << "\t NoIntPointsC:            " << this->es.id.NoIntPointsC << std::endl;
    std::cout << "\t IntPoints:               " << std::endl;
    std::cout << "\t\t\t\t\t\t\t  - "  << "xi:" << "\t\t\t\t" << "eta:" << "\t\t\t" << "zeta:" << "\t\t\t" << "omega:" << std::endl;
    for (int i=0; i<(this->es.id.NoIntPoints); i++)
        std::cout << "\t\t\t\t\t\t\t    "  << std::setw(10) << std::setprecision(9) << std::fixed << this->es.IntPoints[0+4*i] << "   " << this->es.IntPoints[1+4*i] << "   " << this->es.IntPoints[2+4*i] << "   " << this->es.IntPoints[3+4*i] << std::endl;
    return 0;
};
//     int    *NodeSpecs; /* index 0,1,2,.. */
//     double *AdditionalData;
//     char *OutFileName;
//     char *AdditionalNodes;
//     char *AdditionalGraphics;
//     char **CharSwitch;
//     int  *IntSwitch;
//     double *DoubleSwitch;
//     char **SensitivityNames;
//     char *MainTitle;
//     char *SubTitle;
//     char *SubSubTitle;
//     int *SensLowerOrderIndex;/* fi_i, fi_j, fi_k, fi_ij, fi_ik, fi_jk ... for all sens par , all orders */
//     char *MMANextStep;
//     char *MMAStepBack;
//     char *MMAPreIteration;
//     char **IDataNames;
//     int  *IDataIndex;
//     char **RDataNames;
//     int  *RDataIndex;
//     int *ExtraSensitivityData; /*ExtraSensitivityDataLength */
//     char *Bibliography;
//     void *dummy1;
//     void *dummy2;
//     double *PostNodeWeights;
//     int  *MMAInitialisation;
//     double *Version;
//     void *dummy3;
//     void *dummy4;
//     void *dummy5;
//     void *dummy6;
//     void *dummy7;
//     void *dummy8;
//     void *dummy9;
//     void (*Execute)(char *f,char *m);
//     void (*BreakPoint)(char *s, char *c,int m, int p);
//     int (*BoundCheck)(int i,int j,int k,char *m,char *s,char *d,char *e);
// struct{
//         int (*SetElSpec)(struct ElementSpec *,int *,int ,int);
//         void (*SKR)(double *,struct ElementSpec *,ElementData *,NodeSpec **,NodeData **,double *,int *,double *,double **);
//         void (*SRE)(double *,struct ElementSpec *,ElementData *,NodeSpec **,NodeData **,double *,int *,double *);
//         /*forward mode - first order pseudo load*/
//         void (*SSE)(double *,struct ElementSpec *,ElementData *,NodeSpec **,NodeData **,double *,int *,double **);
//         void (*SHI)(double *,struct ElementSpec *,ElementData *,NodeSpec **,NodeData **,double *,int *);
//         void (*SPP)(double *,struct ElementSpec *,ElementData *,NodeSpec **,NodeData **,double *,int *,double **,double **);
//         void (*ELI)(double *,struct ElementSpec *,ElementData *,NodeSpec **,NodeData **,double *,int *);
//         /*position of the nodes at current and previous time step and normal vectors if applicable d[MaxNoNodes][12]*/
//         void (*PAN)(double *,struct ElementSpec *,ElementData *,NodeSpec **,NodeData **,double *,int *,double **);
//         void (*User1)(USER1_ARGUMENTS);
//         void (*User2)(USER2_ARGUMENTS);
//         void (*User3)(USER3_ARGUMENTS);
//         /*backward mode - gradient - first order*/
//         void (*SBG)(double*,struct ElementSpec *, ElementData*, NodeSpec **, NodeData**, double*, int *, void *, double**, double **);
//         /*backward mode - pseudo load*/
//         void (*SBL)(double*,struct ElementSpec *, ElementData*, NodeSpec **, NodeData**, double*, int *, void *, double **);
//         /*forward mode - second order pseudo load*/
//         void(*SS2)(double *, struct ElementSpec *, ElementData *, NodeSpec **, NodeData **, double *, int *, double **);
//         /*backward mode - Hessian - second order*/
//         void (*SBH)(double*,struct ElementSpec *, ElementData*, NodeSpec **, NodeData**, double*, int *, void *, double**, double **);
//         void (*Tasks)(double *,struct ElementSpec *,ElementData *,NodeSpec **,NodeData **,double *,int *,double *,double **,
//               int *,int *,int *,double *,int *, double *);
//         /* static condensation */
//         void (*CondensationData)(struct ElementSpec *,ElementData *, NodeSpec **,NodeData **, double *rdata,
//            int *idata,double **h, double **ha, double **hb,double **dh);
//     }user;