#include <ace_layer.hpp>
#include <iomanip>

int AceGenElement::load(std::string path_to_elememt)
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

int AceGenElement::SKR(
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

int AceGenElement::unload()
{
    if (shared_elmt_lib!=nullptr) dlclose(shared_elmt_lib);  
    if (CMakeDebugConfig) std::cout << "AceGenElement unloaded!" << std::endl;
    return 0;
};

AceGenElement::AceGenElement(std::string path_to_elememt) : AceGenElement()
{
    init(path_to_elememt);
};

AceGenElement::~AceGenElement()
{
    unload();
    delete this->v;
};

int AceGenElement::init(std::string path_to_elememt)
{
    load(path_to_elememt);
    this -> v = new double [this->es.id.WorkingVectorSize];
    return 0;
};

int AceGenElement::status()
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
    std::cout << "\t LocalReKe:               " << this->es.id.LocalReKe << std::endl;
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