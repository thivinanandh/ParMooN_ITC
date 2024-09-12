// =======================================================================
//
// Purpose:     main program for Radiative Transfer Equation equation with new kernels of ParMooN
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 27.04.19
// =======================================================================

#include <Domain.h>
#include <Database.h>
#include <SystemPBE3D.h>
#include <FEDatabase2D.h>
#include <FEDatabase3D.h>
#include <FESpace3D.h>
#include <Output3D.h>
#include <LinAlg.h>

#include <MainUtilities.h>
#include <TimeDiscRout.h>

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <stdlib.h>
#include <MacroCell.h>

// Thivin
#include <SquareStructure1D.h>
#include <SquareMatrix1D.h>
#include <ADISystem1D.h>

// For internal System setup
#include <SystemADI.h>
#include <SystemADI1D.h>

#ifdef _SMPI
#include "mpi.h"
#include <MeshPartition.h>
#endif

// Added for VTK Visualisation
#include <string>
#include <sstream>
#include <iomanip>



double bound = 0;
double timeC = 0;
// =======================================================================
// include current example
// =======================================================================
// #include "../Examples/TCD_3D/rte.h"
// #include "../Examples/TCD_3D/rteSpatial.h"
#include "../Examples/TCD_3D/Siminhale_PBE.h"


void writeVtkFile(const std::string& VtkBaseName,  int img, TOutput3D* Output) {
    std::ostringstream os;
    
    os << "VTK/" << VtkBaseName << "." 
       << std::setfill('0') << std::setw(5) << img << ".vtk";
    Output->WriteVtk(os.str().c_str());
    
}

int main(int argc, char *argv[])
{
    // ======================================================================
    // variable declaration
    // ======================================================================
    int i, j, l, m, N_SubSteps, ORDER, LEVELS, mg_level, N_Cells, N_DOF, img = 1;
    int N_Active, mg_type, N_DOF_Solid, N_Solid_Levels, N_DoF_Intl;

    double *sol, *oldsol, *sol_tmp, *rhs, *oldrhs, t1, t2, errors[5], **Sol_array, **Rhs_array;
    double *Sol_SpatialSurface, *oldsol_all;
    double tau, end_time, *defect, *olderror, *olderror1, *H1T, *L2T, hmin, hmax;
    double start_time, stop_time, start_vtk = 0, end_vtk = 0, total_vtk = 0, l2, H1, L2error_Max = 0., L2error_Max_t;
    double start_assembling = 0, end_assembling = 0, total_assembling = 0;
    double start_solve = 0, end_solve = 0, total_solve = 0;
    double start_int = 0, end_int = 0, total_int = 0;
    double *InternalPts, *Mat_Intl, val, *InternalScaling;

    TDomain *Domain, *Domain_Solid;
    TDatabase *Database = new TDatabase();

    bool UpdateStiffnessMat, UpdateOnlyRhs, ConvectionFirstTime;

    const char vtkdir[] = "VTK";
    char *PsBaseName, *VtkBaseName, *GEO, *PRM, *GEO_SOLID, *PRM_SOLID;
    char Name[] = "name";
    char Description[] = "description";
    char CString[] = "C";
    double Linfty = 0;
    char SubID[] = "";

    int profiling;
#ifdef _SMPI
    int rank, size, out_rank;
    int MaxCpV, MaxSubDomainPerDof;

    MPI_Comm Comm = MPI_COMM_WORLD;
    MPI_Init(&argc, &argv);

    double time1, time2, temp, reduced_errors[4];

    MPI_Comm_rank(Comm, &rank);
    MPI_Comm_size(Comm, &size);

    TDatabase::ParamDB->Comm = Comm;

    int Refine;
    int metisType[2] = {0, 0};
#endif

    TFEDatabase3D *FEDatabase = new TFEDatabase3D();
    TFEDatabase2D *FEDatabase_surface = new TFEDatabase2D();
    TCollection *coll;
    TFESpace3D **Scalar_FeSpaces, *fesp[1], **Scalar_Solid_FeSpaces;
    TFEFunction3D *Scalar_FeFunction, **Scalar_FeFunctions, *Scalar_Solid_FeFunction, **Scalar_Solid_FeFunctions;
    TOutput3D *Output;
    TSystemPBE1D *Internal_System_Space;
    TDomain* Domain_Intl;
    TAuxParam3D *aux;
    MultiIndex3D AllDerivatives[4] = {D000, D100, D010, D001};

    std::ostringstream os;
    os << " ";

    mkdir(vtkdir, 0777);

    // ======================================================================
    // set the database values and generate mesh
    // ======================================================================
    /** set variables' value in TDatabase using argv[1] (*.dat file), and generate the MESH based */
    Domain = new TDomain(argv[1]);

#ifdef _SMPI
    if (rank == 0)
    {
        TDatabase::ParamDB->Par_P0 = 1;
    }
    else
    {
        TDatabase::ParamDB->Par_P0 = 0;
    }
#endif

    OpenFiles();
    OutFile.setf(std::ios::scientific);

#ifdef _SMPI
    if (TDatabase::ParamDB->Par_P0)
#endif
    {
        Database->WriteParamDB(argv[0]);
        Database->WriteTimeDB();
        ExampleFile();
    }

    /* include the mesh from a meshgenerator, for a standard mesh use the build-in function */
    // standard mesh
    PRM = TDatabase::ParamDB->BNDFILE;
    GEO = TDatabase::ParamDB->GEOFILE;
    PsBaseName = TDatabase::ParamDB->BASENAME;
    VtkBaseName = TDatabase::ParamDB->VTKBASENAME;
    /* meshgenerator */
    if (TDatabase::ParamDB->MESH_TYPE == 0)
    {
        Domain->Init(TDatabase::ParamDB->BNDFILE, TDatabase::ParamDB->GEOFILE);
    } // ParMooN  build-in Geo mesh
    else if (TDatabase::ParamDB->MESH_TYPE == 1)
    {
        Domain->GmshGen(TDatabase::ParamDB->GEOFILE);
    } // gmsh mesh
    // else if(TDatabase::ParamDB->MESH_TYPE==2)
    // {
    // Domain->TetrameshGen(TDatabase::ParamDB->GEOFILE);
    // } //tetgen mesh
    else
    {
        OutPut("Mesh Type not known, set MESH_TYPE correctly!!!" << endl);
        exit(0);
    }

    LEVELS = TDatabase::ParamDB->LEVELS;

    if (TDatabase::ParamDB->SOLVER_TYPE == DIRECT)
    {
        TDatabase::ParamDB->UNIFORM_STEPS += (LEVELS - 1);
        LEVELS = 1;
    }
    // refine grid up to the coarsest level
    for (i = 0; i < TDatabase::ParamDB->UNIFORM_STEPS; i++)
        Domain->RegRefineAll();

    if (TDatabase::ParamDB->WRITE_PS)
    {
        // write grid into an Postscript file
        os.seekp(std::ios::beg);
#ifdef _SMPI
        if (TDatabase::ParamDB->Par_P0)
#endif
            os << "Domain" << ".ps" << ends;
        Domain->PS(os.str().c_str(), It_Finest, 0);
    }

    //=========================================================================
    // set data for multigrid
    //=========================================================================
    // set type of multilevel
    mg_type = TDatabase::ParamDB->SC_MG_TYPE_SCALAR;

    if (TDatabase::ParamDB->SOLVER_TYPE == AMG_SOLVE || TDatabase::ParamDB->SOLVER_TYPE == DIRECT)
    {
        mg_type = 0;
        TDatabase::ParamDB->SC_MG_TYPE_SCALAR = mg_type;
    }

    if (mg_type)
    {
        mg_level = LEVELS + 1;
        ORDER = -1;
    }
    else
    {
        mg_level = LEVELS;
        ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
    }

    if (TDatabase::ParamDB->SOLVER_TYPE == GMG)
    {
#ifdef _SMPI
        if (TDatabase::ParamDB->Par_P0)
#endif
        {
            OutPut("=======================================================" << endl);
            OutPut("======           GEOMETRY  LEVEL ");
            OutPut(LEVELS << "              ======" << endl);
            OutPut("======           MULTIGRID LEVEL ");
            OutPut(mg_level << "              ======" << endl);
            OutPut("=======================================================" << endl);
        }
    }

    Scalar_FeSpaces = new TFESpace3D *[mg_level];
    Scalar_FeFunctions = new TFEFunction3D *[mg_level];
    Sol_array = new double *[mg_level];
    Rhs_array = new double *[mg_level];

    //=========================================================================
    // construct all finite element spaces
    // loop over all levels (not a multigrid level but for convergence study)
    //=========================================================================

    for (i = 0; i < LEVELS; i++)
    {
        if (i)
        {
            Domain->RegRefineAll();
        }

        coll = Domain->GetCollection(It_Finest, 0);

        // fespaces for scalar equation
        Scalar_FeSpaces[i] = new TFESpace3D(coll, Name, Description, BoundCondition, ORDER);

        // multilevel multigrid disc
        if (i == LEVELS - 1 && mg_type == 1)
        {
            ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
            Scalar_FeSpaces[mg_level - 1] = new TFESpace3D(coll, Name, Description, BoundCondition, ORDER);
        } //  if(i==LEVELS-1 && i!=mg_level-1)

        //======================================================================
        // construct all finite element functions
        //======================================================================
        N_DOF = Scalar_FeSpaces[i]->GetN_DegreesOfFreedom();
        sol = new double[N_DOF];
        rhs = new double[N_DOF];
        Sol_array[i] = sol;
        Rhs_array[i] = rhs;

        if (i == LEVELS - 1 && mg_type == 1)
        {
            N_DOF = Scalar_FeSpaces[mg_level - 1]->GetN_DegreesOfFreedom();
            sol = new double[N_DOF];
            rhs = new double[N_DOF];
            Sol_array[mg_level - 1] = sol;
            Rhs_array[mg_level - 1] = rhs;

            Scalar_FeFunction = new TFEFunction3D(Scalar_FeSpaces[mg_level - 1], CString, CString, sol, N_DOF);
            Scalar_FeFunctions[mg_level - 1] = Scalar_FeFunction;
        } //   if(i==LEVELS-1 && mg_type==1)
        else
        {
            Scalar_FeFunction = new TFEFunction3D(Scalar_FeSpaces[i], CString, CString, sol, N_DOF);
            Scalar_FeFunctions[i] = Scalar_FeFunction;
        }

#ifdef _MPI
        N_Cells = coll->GetN_Cells();
        printf("rank=%d\t N_Cells   : %d\t Dof all   :%d\n", rank, N_Cells, N_DOF);
#endif
    } // for(i=0;i<LEVELS;i++)

    oldrhs = new double[N_DOF];
    oldsol = new double[N_DOF];

#ifndef _MPI
    N_Cells = coll->GetN_Cells();
#ifdef _SMPI
    if (TDatabase::ParamDB->Par_P0)
#endif
    {
        OutPut("[INFO] N_Cells   : " << N_Cells << endl);
        OutPut("[INFO] Dof all   : " << N_DOF << endl);
        OutPut(endl);
    }
#endif

    //=========================================================================
    // Spatial position at which the internal system is solverd
    //=========================================================================
    double *PosX;
    PosX = new double[3 * N_DOF];
    InternalScaling = new double[N_DOF];
    double *OldIntlRhs = new double[N_DOF];
    double fact;
    TFEVectFunct3D *Pos_VectFunction = new TFEVectFunct3D(Scalar_FeSpaces[mg_level - 1], "Pos", "Pos",
                                                          PosX, N_DOF, 3);
    Pos_VectFunction->GridToData();

    // --- THIVIN ---
    int N_PhySpacePts = N_DOF;
    
    // Store the Boundary Functions in a pointer
    BoundCond1D *BoundConLminLMax[1]; // Function pointer for storing Boundary conditions for the internal co-ordinate
    DoubleFunctND *GrowthAndB_Nuc[1]; // Function pointer for storing the Growth and Nucleation functions for the internal co-ordinate

    
    int N_Cells_Intl; // Numer of cells in the internal co-ordinate
    double StartEnd_intl[2]; // Start and End of the internal co-ordinate
    double** Xpos_internal = new double*[1];

    GetLExampleData(&N_Cells_Intl, StartEnd_intl, BoundConLminLMax, GrowthAndB_Nuc,  Xpos_internal);

    // Create the internal System object
    TSystemADI1D *Internal_System = new TSystemADI1D(N_Cells_Intl, StartEnd_intl[0], StartEnd_intl[1], BoundConLminLMax[0],
                                  GrowthAndB_Nuc[0], Xpos_internal[0]);


    // // Stystem setup for Internal Co-ordinates
    // Internal_System_Space = new TSystemPBE1D(mg_level, Scalar_FeSpaces, Sol_array, Rhs_array,
    //                               TDatabase::ParamDB->DISCTYPE, TDatabase::ParamDB->SOLVER_TYPE, 
    //                               Scalar_FeFunctions, GetKernel, GetRhs, Domain_Intl);  


    int* N_Nodal_position_array = new int[1];
    N_Nodal_position_array[0] =Internal_System->GetN_NodalPts();
    double **nodal_positions = new double*[1];
    nodal_positions[0] = Internal_System->GetNodalPt_Coord();

    int N_Nodals_All = N_PhySpacePts * N_Nodal_position_array[0];
    double* solution_all = new double[N_Nodals_All];

    int N_Coord = 3 ;// (x,y,z )
    Internal_System->Init(N_Coord, N_PhySpacePts, PosX, 1, N_Nodal_position_array, nodal_positions, 0, solution_all);

    int N_InternalPts = N_Nodal_position_array[0];

    // print the Stats of the internal system
    cout << "Internal System Stats: " << endl;
    cout << "Num Internal Nodal Points: " << N_Nodal_position_array[0] << endl;
    cout << "Num Physical Space Points: " << N_PhySpacePts << endl;
    cout << "Num Total Nodal Points: " << N_Nodals_All << endl;


    // Interpolate the Initial Solution
    Internal_System->Interpolate_With_Coord(N_Coord, solution_all, InitialValue);

    // -- OUTPUT SETUP FOR VISUALISATION ---- //

    // Create an double array to store the solution of the internal system
    double** solution_visualize = new double*[N_InternalPts];
    for (int i = 0; i < N_InternalPts; ++i) {
        solution_visualize[i] = new double[N_PhySpacePts]();
    }

    // Generate Individual FEFunctions for the internal system
    TFEFunction3D **Scalar_FeFunctions_Intl = new TFEFunction3D*[N_InternalPts];
    for (int i = 0; i < N_InternalPts; i++)
    {
        char name[10];
        sprintf(name, "Intl_%d", i);
        Scalar_FeFunctions_Intl[i] = new TFEFunction3D(Scalar_FeSpaces[mg_level - 1], name, name, solution_visualize[i], N_PhySpacePts);
    }
    


    // ---------- Obtain the alternate order of the solution for VTK output 
    // Lets say there are 9 Physical Points (3D) (Identified as A, B, C ...) and 2 internal points(identified as l1, l2) per physical point
    // The Solution in a 3D + 1D is represented as Al1, Al2, Bl1, Bl2, Cl1, Cl2, Dl1, Dl2, El1, El2, Fl1, Fl2, Gl1, Gl2, Hl1, Hl2, Il1, Il2
    // The current solution is saved as 9*2 = 18, where the array is stored as Al1, Al2, Bl1, Bl2, Cl1, Cl2, Dl1, Dl2, El1, El2, Fl1, Fl2, Gl1, Gl2, Hl1, Hl2, Il1, Il2


    // However, for visualization, we need to store the solution as Al1, Bl1, Cl1, Dl1, El1, Fl1, Gl1, Hl1, Il1, Al2, Bl2, Cl2, Dl2, El2, Fl2, Gl2, Hl2, Il2
    // This is done by the following code
    
    // // copy the corresponding levels solution to the solution_visualize array
    for (int i = 0; i < N_PhySpacePts * N_InternalPts; ++i) {
        int old_index = (i % N_PhySpacePts) * N_InternalPts + (i / N_PhySpacePts);
        int level = i / N_PhySpacePts;
        solution_visualize[level][i % N_PhySpacePts] = solution_all[old_index];
    }
    
    
    // Output for Visualization of the codes. 
    VtkBaseName = TDatabase::ParamDB->VTKBASENAME;

    TOutput3D* Output_Intl = new TOutput3D(2, 2, 1, 1, Domain);
    
    // Add all the internal levels for visualisation
    for (int i = 0; i < N_InternalPts; i++)
        Output_Intl->AddFEFunction(Scalar_FeFunctions_Intl[i]);

    
    // Output the solution for all internal layers
    writeVtkFile(VtkBaseName, i, Output_Intl);
    img++;
    //======================================================================
    // parameters for time stepping scheme
    //======================================================================  
    m = 0;
    N_SubSteps = GetN_SubSteps();
    end_time = TDatabase::TimeDB->ENDTIME;

    // Assemble all the Internal System matrices
    Internal_System->AssembleMassMat();
}