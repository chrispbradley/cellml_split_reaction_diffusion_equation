PROGRAM CellMLSplitReactionDiffusionEquation

  USE OpenCMISS
  USE OpenCMISS_Iron
#ifndef NOMPIMOD
  USE MPI
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif

  !-----------------------------------------------------------------------------------------------------------
  ! PROGRAM VARIABLES AND TYPES
  !-----------------------------------------------------------------------------------------------------------

  !Test program parameters

  REAL(CMISSRP), PARAMETER :: LENGTH=100.0_CMISSRP
  
  INTEGER(CMISSIntg), PARAMETER :: CONTEXT_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: COORDINATE_SYSTEM_USER_NUMBER=2
  INTEGER(CMISSIntg), PARAMETER :: REGION_USER_NUMBER=3
  INTEGER(CMISSIntg), PARAMETER :: BASIS_USER_NUMBER=4
  INTEGER(CMISSIntg), PARAMETER :: GENERATED_MESH_USER_NUMBER=5
  INTEGER(CMISSIntg), PARAMETER :: MESH_USER_NUMBER=6
  INTEGER(CMISSIntg), PARAMETER :: DECOMPOSITION_USER_NUMBER=7
  INTEGER(CMISSIntg), PARAMETER :: DECOMPOSER_USER_NUMBER=8
  INTEGER(CMISSIntg), PARAMETER :: GEOMETRIC_FIELD_USER_NUMBER=9
  INTEGER(CMISSIntg), PARAMETER :: EQUATIONS_SET_FIELD_USER_NUMBER=10
  INTEGER(CMISSIntg), PARAMETER :: DEPENDENT_FIELD_USER_NUMBER=11
  INTEGER(CMISSIntg), PARAMETER :: MATERIALS_FIELD_USER_NUMBER=12
  INTEGER(CMISSIntg), PARAMETER :: SOURCE_FIELD_USER_NUMBER=13
  INTEGER(CMISSIntg), PARAMETER :: CELLML_MODELS_FIELD_USER_NUMBER=14
  INTEGER(CMISSIntg), PARAMETER :: CELLML_STATE_FIELD_USER_NUMBER=15
  INTEGER(CMISSIntg), PARAMETER :: CELLML_INTERMEDIATE_FIELD_USER_NUMBER=16
  INTEGER(CMISSIntg), PARAMETER :: CELLML_PARAMETERS_FIELD_USER_NUMBER=17
  INTEGER(CMISSIntg), PARAMETER :: EQUATIONS_SET_USER_NUMBER=18
  INTEGER(CMISSIntg), PARAMETER :: CELLML_USER_NUMBER=19
  INTEGER(CMISSIntg), PARAMETER :: PROBLEM_USER_NUMBER=20

  !Program types
  
  !Program variables

  INTEGER(CMISSIntg) :: numberOfGlobalXElements
  INTEGER(CMISSIntg) :: nodeDomain,nodeIdx,nodeNumber
  INTEGER(CMISSIntg) :: boundaryConditionNodes(2)  
  INTEGER(CMISSIntg) :: constantModelIndex
  
  !CMISS variables

  TYPE(cmfe_BasisType) :: basis
  TYPE(cmfe_BoundaryConditionsType) :: boundaryConditions
  TYPE(cmfe_CellMLType) :: cellML
  TYPE(cmfe_CellMLEquationsType) :: cellMLEquations
  TYPE(cmfe_ComputationEnvironmentType) :: computationEnvironment
  TYPE(cmfe_ContextType) :: context
  TYPE(cmfe_ControlLoopType) :: controlLoop
  TYPE(cmfe_CoordinateSystemType) :: coordinateSystem
  TYPE(cmfe_DecompositionType) :: decomposition
  TYPE(cmfe_DecomposerType) :: decomposer
  TYPE(cmfe_EquationsType) :: equations
  TYPE(cmfe_EquationsSetType) :: equationsSet
  TYPE(cmfe_FieldType) :: geometricField,equationsSetField,dependentField,materialsField,sourceField
  TYPE(cmfe_FieldType) :: cellMLModelsField,cellMLStateField,cellMLIntermediateField,cellMLParametersField
  TYPE(cmfe_FieldsType) :: fields
  TYPE(cmfe_GeneratedMeshType) :: generatedMesh  
  TYPE(cmfe_MeshType) :: mesh
  TYPE(cmfe_ProblemType) :: problem
  TYPE(cmfe_RegionType) :: region,worldRegion
  TYPE(cmfe_SolverType) :: solver, linearSolver
  TYPE(cmfe_SolverEquationsType) :: solverEquations
  TYPE(cmfe_WorkGroupType) :: worldWorkGroup

  LOGICAL :: exportField
  
  !Generic CMISS variables
  
  INTEGER(CMISSIntg) :: numberOfComputationalNodes,computationalNodeNumber
  INTEGER(CMISSIntg) :: decompositionIndex,equationsSetIndex,cellMLIndex
  INTEGER(CMISSIntg) :: err

  !-----------------------------------------------------------------------------------------------------------
  ! PROBLEM CONTROL PANEL
  !-----------------------------------------------------------------------------------------------------------

  !Intialise OpenCMISS
  CALL cmfe_Initialise(err)
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,err)
  !Create a context
  CALL cmfe_Context_Initialise(context,err)
  CALL cmfe_Context_Create(CONTEXT_USER_NUMBER,context,err)
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,err)
  CALL cmfe_Region_Initialise(worldRegion,err)
  CALL cmfe_Context_WorldRegionGet(context,worldRegion,err)
  
  !Get the computational nodes information
  CALL cmfe_ComputationEnvironment_Initialise(computationEnvironment,err)
  CALL cmfe_Context_ComputationEnvironmentGet(context,computationEnvironment,err)
  
  CALL cmfe_WorkGroup_Initialise(worldWorkGroup,err)
  CALL cmfe_ComputationEnvironment_WorldWorkGroupGet(computationEnvironment,worldWorkGroup,err)
  CALL cmfe_WorkGroup_NumberOfGroupNodesGet(worldWorkGroup,numberOfComputationalNodes,err)
  CALL cmfe_WorkGroup_GroupNodeNumberGet(worldWorkGroup,computationalNodeNumber,err)

  numberOfGlobalXElements=10

  !-----------------------------------------------------------------------------------------------------------
  ! COORDINATE SYSTEM
  !-----------------------------------------------------------------------------------------------------------  

  !Start the creation of a new RC coordinate system
  CALL cmfe_CoordinateSystem_Initialise(coordinateSystem,err)
  CALL cmfe_CoordinateSystem_CreateStart(COORDINATE_SYSTEM_USER_NUMBER,context,coordinateSystem,err)
  !Set the coordinate system to be 1D
  CALL cmfe_CoordinateSystem_DimensionSet(coordinateSystem,1,err)
  !Finish the creation of the coordinate system
  CALL cmfe_CoordinateSystem_CreateFinish(coordinateSystem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! REGION
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of the region
  CALL cmfe_Region_Initialise(region,err)
  CALL cmfe_Region_CreateStart(REGION_USER_NUMBER,worldRegion,region,err)
  CALL cmfe_Region_LabelSet(region,"cellml_split_reaction_diffusion_equation",err)
  !Set the regions coordinate system to the 1D RC coordinate system that we have created
  CALL cmfe_Region_CoordinateSystemSet(region,coordinateSystem,err)
  !Finish the creation of the region
  CALL cmfe_Region_CreateFinish(region,err)

  !-----------------------------------------------------------------------------------------------------------
  ! BASIS
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of a basis (default is linear lagrange)
  CALL cmfe_Basis_Initialise(basis,err)
  CALL cmfe_Basis_CreateStart(BASIS_USER_NUMBER,context,basis,err)
  !Set the basis to be a linear Lagrange basis
  CALL cmfe_Basis_NumberOfXiSet(basis,1,err)
  !Finish the creation of the basis
  CALL cmfe_Basis_CreateFinish(basis,err)

  !-----------------------------------------------------------------------------------------------------------
  ! MESH
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of a generated mesh in the region
  CALL cmfe_GeneratedMesh_Initialise(generatedMesh,err)
  CALL cmfe_GeneratedMesh_CreateStart(GENERATED_MESH_USER_NUMBER,region,generatedMesh,err)
  !Set up a regular 1D mesh
  CALL cmfe_GeneratedMesh_TypeSet(generatedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,err)
  !Set the default basis
  CALL cmfe_GeneratedMesh_BasisSet(generatedMesh,basis,err)   
  !Define the mesh on the region
  CALL cmfe_GeneratedMesh_ExtentSet(generatedMesh,[LENGTH],err)
  CALL cmfe_GeneratedMesh_NumberOfElementsSet(generatedMesh,[numberOfGlobalXElements],err)
  !Finish the creation of a generated mesh in the region
  CALL cmfe_Mesh_Initialise(mesh,err)
  CALL cmfe_GeneratedMesh_CreateFinish(generatedMesh,MESH_USER_NUMBER,mesh,err)

  !-----------------------------------------------------------------------------------------------------------
  ! DECOMPOSITION
  !-----------------------------------------------------------------------------------------------------------

  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(decomposition,err)
  CALL cmfe_Decomposition_CreateStart(DECOMPOSITION_USER_NUMBER,mesh,decomposition,err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL cmfe_Decomposition_TypeSet(decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,err)
  !Finish the decomposition
  CALL cmfe_Decomposition_CreateFinish(decomposition,err)

  !-----------------------------------------------------------------------------------------------------------
  ! DECOMPOSER
  !-----------------------------------------------------------------------------------------------------------

  CALL cmfe_Decomposer_Initialise(decomposer,err)
  CALL cmfe_Decomposer_CreateStart(DECOMPOSER_USER_NUMBER,region,worldWorkGroup,decomposer,err)
  !Add in the decomposition
  CALL cmfe_Decomposer_DecompositionAdd(decomposer,decomposition,decompositionIndex,err)
  !Finish the decomposer
  CALL cmfe_Decomposer_CreateFinish(decomposer,err)
  
  !-----------------------------------------------------------------------------------------------------------
  ! GEOMETRIC FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Start to create a default (geometric) field on the region
  CALL cmfe_Field_Initialise(geometricField,err)
  CALL cmfe_Field_CreateStart(GEOMETRIC_FIELD_USER_NUMBER,region,geometricField,err)
  !Set the decomposition to use
  CALL cmfe_Field_DecompositionSet(geometricField,decomposition,err)
  !Set the domain to be used by the field components.
  CALL cmfe_Field_ComponentMeshComponentSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,err)
  !Finish creating the field
  CALL cmfe_Field_CreateFinish(geometricField,err)

  !Update the geometric field parameters
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(generatedMesh,geometricField,err)

  !-----------------------------------------------------------------------------------------------------------
  ! EQUATIONS SETS
  !-----------------------------------------------------------------------------------------------------------

  !Create the cellml reaction with split reaction diffusion equations_set 
  CALL cmfe_EquationsSet_Initialise(equationsSet,err)
  CALL cmfe_Field_Initialise(equationsSetField,err)
  CALL cmfe_EquationsSet_CreateStart(EQUATIONS_SET_USER_NUMBER,region,geometricField,[CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMFE_EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE,CMFE_EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE], &
    & EQUATIONS_SET_FIELD_USER_NUMBER,equationsSetField,equationsSet,err)
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(equationsSet,err)

  !-----------------------------------------------------------------------------------------------------------
  ! DEPENDENT FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set dependent field variables (primary and secondary variable)
  CALL cmfe_Field_Initialise(dependentField,err)
  CALL cmfe_EquationsSet_DependentCreateStart(equationsSet,DEPENDENT_FIELD_USER_NUMBER,dependentField,err)
  CALL cmfe_Field_VariableLabelSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,"U",err)
  CALL cmfe_Field_VariableLabelSet(dependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,"DELUDELN",err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(equationsSet,err)

  !-----------------------------------------------------------------------------------------------------------
  ! MATERIAL FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set material field variables
  !by default 2 components for reaction-diffusion i.e. diffusion coefficient in x direction 
  !set constant spatially = 1, and storage coefficient set to 1
  CALL cmfe_Field_Initialise(materialsField,err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(equationsSet,MATERIALS_FIELD_USER_NUMBER,materialsField,err)
  CALL cmfe_Field_VariableLabelSet(materialsField,CMFE_FIELD_U_VARIABLE_TYPE,"Material",err)
  !Finish the equations set materials field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(equationsSet,err)
  CALL cmfe_Field_ComponentValuesInitialise(materialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
   & 1,0.5_CMISSRP,err) !diffusion coefficent in x
  CALL cmfe_Field_ComponentValuesInitialise(materialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
   & 2,1.0_CMISSRP,err) ! storage coefficient
  
  !-----------------------------------------------------------------------------------------------------------
  ! SOURCE FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Set up source field for reaction diffusion equation set. Note that for the split problem subtype, 
  !the source field is not used at all.
  CALL cmfe_Field_Initialise(sourceField,err)
  CALL cmfe_EquationsSet_SourceCreateStart(equationsSet,SOURCE_FIELD_USER_NUMBER,sourceField,err)
  CALL cmfe_Field_VariableLabelSet(sourceField,CMFE_FIELD_U_VARIABLE_TYPE,"Source",err)
  !Finish the equations set source field variables
  CALL cmfe_EquationsSet_SourceCreateFinish(equationsSet,err)
  CALL cmfe_Field_ComponentValuesInitialise(sourceField,CMFE_FIELD_U_VARIABLE_TYPE, &
    & CMFE_FIELD_VALUES_SET_TYPE,1,0.0_CMISSRP,err)

  !-----------------------------------------------------------------------------------------------------------
  ! CELLML FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Start to set up CellML Fields

  !Create the CellML environment
  CALL cmfe_CellML_Initialise(cellML,err)
  CALL cmfe_CellML_CreateStart(CELLML_USER_NUMBER,region,cellML,err)
  !Import a constant source (i.e. rate of generation/depletion is zero) model from a file
  CALL cmfe_CellML_ModelImport(cellML,"constant_rate.xml",constantModelIndex,err)
  !Speify the variables in the imported model that will be used. 
  
  ! Now we have imported all the models we are able to specify which variables from the model we want:
  ! - to set from this side
  !These are effectively parameters that won't change in the course of the ode solving for one time step. 
  !i.e. fixed before running cellml, known in opencmiss and changed only in opencmiss - components of the parameters field
  CALL cmfe_CellML_VariableSetAsKnown(cellML,constantModelIndex,"dude/param",err)
  ! - to get from the CellML side. variables in cellml model that are not state variables, but are dependent on 
  !independent and state variables. - components of intermediate field
  CALL cmfe_CellML_VariableSetAsWanted(cellML,constantModelIndex,"dude/intmd",err)
  !Finish the CellML environment
  CALL cmfe_CellML_CreateFinish(cellML,err)

  !Start the creation of CellML <--> OpenCMISS field maps
  CALL cmfe_CellML_FieldMapsCreateStart(cellML,err)
  !Set up the field variable component <--> CellML model variable mappings.
  !Here opencmiss fields are mapped to appropriate cellml fields (parameters/intermediates/state).
  !in monodomain problems, typically one wants Vm, a state variable. In the present case, Ca concentration is the 
  !state variable. Furthermore, an injection of Ca is also required. Since in this example 'order-splitting' is used, 
  !the source field is not required. Instead, the problem is solved as an ODE and then as a PDE separately.
  !Ca concentration is obtained from the dependent field and solve the DAE to give new values of the Ca concentration
  !in the dependent field again. This is then used as the initial set up for the dynamic solver. Thus, the dependent field 
  !is mapped to appropriate component variable names.cellml/opencmiss will look up the appropriate field.

  !On the other hand, if there is no order-splitting, i.e. a proper reaction-diffusion equation to solve, then the current Ca
  !concentration from the dependent field is used to solve the DAE and the result of the DAE is substituted back to 
  !the source field.
 
  CALL cmfe_CellML_CreateFieldToCellMLMap(cellML,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
    & constantModelIndex,"dude/ca",CMFE_FIELD_VALUES_SET_TYPE,err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(cellML,constantModelIndex,"dude/ca",CMFE_FIELD_VALUES_SET_TYPE, &
    & dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,err)
  !Finish the creation of CellML <--> OpenCMISS field maps
  CALL cmfe_CellML_FieldMapsCreateFinish(cellML,err)

  !set initial value of the dependent field/state variable, Ca concentration.
  CALL cmfe_Field_ComponentValuesInitialise(dependentField,CMFE_FIELD_U_VARIABLE_TYPE, &
    & CMFE_FIELD_VALUES_SET_TYPE,1,0.0_CMISSRP,err)
  nodeIdx=2
  CALL cmfe_Decomposition_NodeDomainGet(decomposition,nodeIdx,1,nodeDomain,err)
  IF(nodeDomain==computationalNodeNumber) THEN
    CALL cmfe_Field_ParameterSetUpdateNode(dependentField,CMFE_FIELD_U_VARIABLE_TYPE, &
     & CMFE_FIELD_VALUES_SET_TYPE, &
     & 1,1,nodeIdx,1,0.0_CMISSRP,err) 
  ENDIF
  !Start the creation of the CellML models field. This field is an integer field that stores which nodes have which cellml 
  !model
  CALL cmfe_Field_Initialise(cellMLModelsField,err)
  CALL cmfe_CellML_ModelsFieldCreateStart(cellML, CELLML_MODELS_FIELD_USER_NUMBER, &
    & cellMLModelsField,err)
  !Finish the creation of the CellML models field
  CALL cmfe_CellML_ModelsFieldCreateFinish(cellML,err)
  !The CellMLModelsField is an integer field that stores which model is being used by which node.
  !By default all field parameters have default model value of 1, i.e. the first model. But, this command below is 
  !for example purposes
  CALL cmfe_Field_ComponentValuesInitialise(cellMLModelsField,CMFE_FIELD_U_VARIABLE_TYPE, &
    & CMFE_FIELD_VALUES_SET_TYPE,1,1_CMISSIntg,err)

  !Set up the models field
  !DO N=1,(numberOfGlobalXElements+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)*(NUMBER_GLOBAL_Z_ELEMENTS+1)
  !  IF(N < 5) THEN
  !    CELL_TYPE = 1
  !  ELSE
  !    CELL_TYPE = 2
  !  ENDIF
  !  CALL cmfe_FieldParameterSetUpdateNode(cellMLModelsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,N,1,CELL_TYPE,err)
  !END DO
  !CALL cmfe_FieldParameterSetUpdateStart(cellMLModelsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,err)
  !CALL cmfe_FieldParameterSetUpdateFinish(cellMLModelsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,err)

  !Start the creation of the cellML state field
  CALL cmfe_Field_Initialise(cellMLStateField,err)
  CALL cmfe_CellML_StateFieldCreateStart(cellML,CELLML_STATE_FIELD_USER_NUMBER,cellMLStateField,err)
  !Finish the creation of the CellML state field
  CALL cmfe_CellML_StateFieldCreateFinish(cellML,err)

  !Start the creation of the CellML intermediate field
  CALL cmfe_Field_Initialise(cellMLIntermediateField,err)
  CALL cmfe_CellML_IntermediateFieldCreateStart(cellML,CELLML_INTERMEDIATE_FIELD_USER_NUMBER,cellMLIntermediateField,err)
  !Finish the creation of the CellML intermediate field
  CALL cmfe_CellML_IntermediateFieldCreateFinish(cellML,err)

  !Start the creation of CellML parameters field
  CALL cmfe_Field_Initialise(cellMLParametersField,err)
  CALL cmfe_CellML_ParametersFieldCreateStart(cellML,CELLML_PARAMETERS_FIELD_USER_NUMBER,cellMLParametersField,err)
  !Finish the creation of CellML parameters
  CALL cmfe_CellML_ParametersFieldCreateFinish(cellML,err)

  !-----------------------------------------------------------------------------------------------------------
  ! EQUATIONS SETS
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set equations
  CALL cmfe_Equations_Initialise(equations,err)
  CALL cmfe_EquationsSet_EquationsCreateStart(equationsSet,equations,err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(equations,CMFE_EQUATIONS_SPARSE_MATRICES,err)
  !Set the equations set output
  CALL cmfe_Equations_OutputTypeSet(equations,CMFE_EQUATIONS_NO_OUTPUT,err)
  !CALL cmfe_EquationsOutputTypeSet(equations,CMFE_EQUATIONS_TIMING_OUTPUT,err)
  !CALL cmfe_EquationsOutputTypeSet(equations,CMFE_EQUATIONS_MATRIX_OUTPUT,err)
  !CALL cmfe_equationsOutputTypeSet(equations,cmfe_EQUATIONS_ELEMENT_MATRIX_OUTPUT,err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(equationsSet,err)

  !-----------------------------------------------------------------------------------------------------------
  ! PROBLEM
  !-----------------------------------------------------------------------------------------------------------

  !Create the problem
  CALL cmfe_Problem_Initialise(problem,err)
  CALL cmfe_Problem_CreateStart(PROBLEM_USER_NUMBER,context,[CMFE_PROBLEM_CLASSICAL_FIELD_CLASS, &
    & CMFE_PROBLEM_REACTION_DIFFUSION_EQUATION_TYPE,CMFE_PROBLEM_CELLML_REAC_INTEG_REAC_DIFF_STRANG_SPLIT_SUBTYPE],problem,err)
  !Finish the creation of a problem.
  CALL cmfe_Problem_CreateFinish(problem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! PROBLEM CONTROL
  !-----------------------------------------------------------------------------------------------------------

  !Create the problem control
  CALL cmfe_Problem_ControlLoopCreateStart(problem,err)
  !Get the control loop
  CALL cmfe_ControlLoop_Initialise(controlLoop,err)
  CALL cmfe_Problem_ControlLoopGet(problem,CMFE_CONTROL_LOOP_NODE,controlLoop,err)
  !Set the times
  CALL cmfe_ControlLoop_TimesSet(controlLoop,0.0_CMISSRP,0.5_CMISSRP,0.01_CMISSRP,err)
  CALL cmfe_ControlLoop_OutputTypeSet(controlLoop,CMFE_CONTROL_LOOP_PROGRESS_OUTPUT,err)
  !Finish creating the problem control loop
  CALL cmfe_Problem_ControlLoopCreateFinish(problem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! PROBLEM SOLVERS
  !-----------------------------------------------------------------------------------------------------------

  !Set up the problem solvers for Strang splitting
  CALL cmfe_Problem_SolversCreateStart(problem,err)
  !First solver is a DAE solver
  CALL cmfe_Solver_Initialise(solver,err)
  CALL cmfe_Problem_SolverGet(problem,CMFE_CONTROL_LOOP_NODE,1,solver,err)
  CALL cmfe_Solver_DAESolverTypeSet(solver,CMFE_SOLVER_DAE_EULER,err)
  CALL cmfe_Solver_DAETimeStepSet(solver,0.0000001_CMISSRP,err)
  CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_MATRIX_OUTPUT,err)

  !Second solver is the dynamic solver for solving the parabolic equation
  CALL cmfe_Solver_Initialise(solver,err)
  CALL cmfe_Solver_Initialise(linearSolver,err)
  CALL cmfe_Problem_SolverGet(problem,CMFE_CONTROL_LOOP_NODE,2,solver,err)
  !set theta - backward vs forward time step parameter
  CALL cmfe_Solver_DynamicThetaSet(solver,1.0_CMISSRP,err)
  !CALL cmfe_SolverOutputTypeSet(solver,CMFE_SOLVER_NO_OUTPUT,err)
  CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_TIMING_OUTPUT,err)
  !get the dynamic linear solver from the solver
  CALL cmfe_Solver_DynamicLinearSolverGet(solver,linearSolver,err)
  CALL cmfe_Solver_LibraryTypeSet(linearSolver,CMFE_SOLVER_LAPACK_LIBRARY,err)
  CALL cmfe_Solver_LinearTypeSet(linearSolver,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,err)
  CALL cmfe_Solver_LinearDirectTypeSet(linearSolver,CMFE_SOLVER_DIRECT_LU,err)
  !CALL cmfe_SolverLibraryTypeSet(linearSolver,CMFE_SOLVER_CMISS_LIBRARY,err)
  !CALL cmfe_SolverLinearTypeSet(linearSolver,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,err)
  !CALL cmfe_SolverLibraryTypeSet(linearSolver,CMFE_SOLVER_MUMPS_LIBRARY,err)
  !CALL cmfe_Solver_LinearIterativeMaximumIterationsSet(linearSolver,10000,err)

  !Third solver is another DAE solver
  CALL cmfe_Solver_Initialise(solver,err)
  CALL cmfe_Problem_SolverGet(problem,CMFE_CONTROL_LOOP_NODE,3,solver,err)
  CALL cmfe_Solver_DAESolverTypeSet(solver,CMFE_SOLVER_DAE_EULER,err)
  CALL cmfe_Solver_DAETimeStepSet(solver,0.0000001_CMISSRP,err)
  CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_MATRIX_OUTPUT,err)
  !CALL cmfe_SolverOutputTypeSet(solver,CMFE_SOLVER_TIMING_OUTPUT,err)
  !CALL cmfe_SolverOutputTypeSet(solver,CMFE_SOLVER_SOLVER_OUTPUT,err)
  !CALL cmfe_SolverOutputTypeSet(solver,CMFE_SOLVER_PROGRESS_OUTPUT,err)
  !Finish the creation of the problem solver
  CALL cmfe_Problem_SolversCreateFinish(problem,err)

  !Start the creation of the problem solver CellML equations
  CALL cmfe_Problem_CellMLEquationsCreateStart(problem,err)
  !Get the first solver  
  !Get the CellML equations
  CALL cmfe_Solver_Initialise(solver,err)
  CALL cmfe_Problem_SolverGet(problem,CMFE_CONTROL_LOOP_NODE,1,solver,err)
  CALL cmfe_CellMLEquations_Initialise(cellMLEquations,err)
  CALL cmfe_Solver_CellMLEquationsGet(solver,cellMLEquations,err)
  !Add in the CellML environement
  CALL cmfe_CellMLEquations_CellMLAdd(cellMLEquations,cellML,cellMLIndex,err)

  !Get the third solver  
  !Get the CellML equations
  CALL cmfe_Solver_Initialise(solver,err)
  CALL cmfe_Problem_SolverGet(problem,CMFE_CONTROL_LOOP_NODE,3,solver,err)
  CALL cmfe_CellMLEquations_Initialise(cellMLEquations,err)
  CALL cmfe_Solver_CellMLEquationsGet(solver,cellMLEquations,err)
  !Add in the CellML environement
  CALL cmfe_CellMLEquations_CellMLAdd(cellMLEquations,cellML,cellMLIndex,err)

  !Finish the creation of the problem solver CellML equations
  CALL cmfe_Problem_CellMLEquationsCreateFinish(problem,err)

  !Start the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateStart(problem,err)
  !Get the second solver  
  !Get the solver equations
  CALL cmfe_Solver_Initialise(solver,err)
  CALL cmfe_Problem_SolverGet(problem,CMFE_CONTROL_LOOP_NODE,2,solver,err)
  CALL cmfe_SolverEquations_Initialise(solverEquations,err)
  CALL cmfe_Solver_SolverEquationsGet(solver,solverEquations,err)
  !Set the solver equations sparsity
  !CALL cmfe_SolverEquationsSparsityTypeSet(solverEquations,CMFE_SOLVER_EQUATIONS_SPARSE_MATRICES,err)
  CALL cmfe_SolverEquations_SparsityTypeSet(solverEquations,CMFE_SOLVER_SPARSE_MATRICES,err)  
  !Add in the equations set
  CALL cmfe_SolverEquations_EquationsSetAdd(solverEquations,equationsSet,equationsSetIndex,err)
  !Finish the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateFinish(problem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! BOUNDARY CONDITIONS
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set boundary conditions
  boundaryConditionNodes = [1,11]
  CALL cmfe_BoundaryConditions_Initialise(boundaryConditions,err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(solverEquations,boundaryConditions,err)
  DO nodeIdx=1,2
    nodeNumber = boundaryConditionNodes(nodeIdx)
    CALL cmfe_Decomposition_NodeDomainGet(decomposition,nodeNumber,1,nodeDomain,err)
    IF(nodeDomain==computationalNodeNumber) THEN
      CALL cmfe_BoundaryConditions_SetNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_NO_GLOBAL_DERIV, &
        & nodeNumber,1,CMFE_BOUNDARY_CONDITION_FIXED,1.5_CMISSRP,err)
      !Need to set CellML model to zero at the nodes at which value of ca has been fixed.
      CALL cmfe_Field_ParameterSetUpdateNode(cellMLModelsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,&
        & 1,1,nodeNumber,1,0_CMISSIntg,err) 
    ENDIF
  ENDDO
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(solverEquations,err)

  !-----------------------------------------------------------------------------------------------------------
  ! SOLVER
  !-----------------------------------------------------------------------------------------------------------

  !Solve the problem
  CALL cmfe_Problem_Solve(problem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! OUTPUT
  !-----------------------------------------------------------------------------------------------------------
 
  exportField=.TRUE.
  IF(exportField) THEN
    CALL cmfe_Fields_Initialise(fields,err)
    CALL cmfe_Fields_Create(region,fields,err)
    CALL cmfe_Fields_NodesExport(fields,"cellml_split_reaction_diffusion_equation","FORTRAN",err)
    CALL cmfe_Fields_ElementsExport(fields,"cellml_split_reaction_diffusion_equation","FORTRAN",err)
    CALL cmfe_Fields_Finalise(fields,err)
  ENDIF

  !Destroy the context
  CALL cmfe_Context_Destroy(context,err)
  !Finalise OpenCMISS
  CALL cmfe_Finalise(err)
  
  WRITE(*,'(A)') "Program successfully completed."

  STOP
  
END PROGRAM CellMLSplitReactionDiffusionEquation
