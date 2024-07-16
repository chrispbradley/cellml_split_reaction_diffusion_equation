PROGRAM CellMLSplitReactionDiffusionEquation

  USE OpenCMISS

  IMPLICIT NONE

  !-----------------------------------------------------------------------------------------------------------
  ! PROGRAM VARIABLES AND TYPES
  !-----------------------------------------------------------------------------------------------------------

  !Test program parameters

  REAL(OC_RP), PARAMETER :: LENGTH=100.0_OC_RP
  
  INTEGER(OC_Intg), PARAMETER :: CONTEXT_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: COORDINATE_SYSTEM_USER_NUMBER=2
  INTEGER(OC_Intg), PARAMETER :: REGION_USER_NUMBER=3
  INTEGER(OC_Intg), PARAMETER :: BASIS_USER_NUMBER=4
  INTEGER(OC_Intg), PARAMETER :: GENERATED_MESH_USER_NUMBER=5
  INTEGER(OC_Intg), PARAMETER :: MESH_USER_NUMBER=6
  INTEGER(OC_Intg), PARAMETER :: DECOMPOSITION_USER_NUMBER=7
  INTEGER(OC_Intg), PARAMETER :: DECOMPOSER_USER_NUMBER=8
  INTEGER(OC_Intg), PARAMETER :: GEOMETRIC_FIELD_USER_NUMBER=9
  INTEGER(OC_Intg), PARAMETER :: EQUATIONS_SET_FIELD_USER_NUMBER=10
  INTEGER(OC_Intg), PARAMETER :: DEPENDENT_FIELD_USER_NUMBER=11
  INTEGER(OC_Intg), PARAMETER :: MATERIALS_FIELD_USER_NUMBER=12
  INTEGER(OC_Intg), PARAMETER :: SOURCE_FIELD_USER_NUMBER=13
  INTEGER(OC_Intg), PARAMETER :: CELLML_MODELS_FIELD_USER_NUMBER=14
  INTEGER(OC_Intg), PARAMETER :: CELLML_STATE_FIELD_USER_NUMBER=15
  INTEGER(OC_Intg), PARAMETER :: CELLML_INTERMEDIATE_FIELD_USER_NUMBER=16
  INTEGER(OC_Intg), PARAMETER :: CELLML_PARAMETERS_FIELD_USER_NUMBER=17
  INTEGER(OC_Intg), PARAMETER :: EQUATIONS_SET_USER_NUMBER=18
  INTEGER(OC_Intg), PARAMETER :: CELLML_USER_NUMBER=19
  INTEGER(OC_Intg), PARAMETER :: PROBLEM_USER_NUMBER=20

  !Program types
  
  !Program variables

  INTEGER(OC_Intg) :: numberOfGlobalXElements
  INTEGER(OC_Intg) :: nodeDomain,nodeIdx,nodeNumber
  INTEGER(OC_Intg) :: boundaryConditionNodes(2)  
  INTEGER(OC_Intg) :: constantModelIndex
  
  !CMISS variables

  TYPE(OC_BasisType) :: basis
  TYPE(OC_BoundaryConditionsType) :: boundaryConditions
  TYPE(OC_CellMLType) :: cellML
  TYPE(OC_CellMLEquationsType) :: cellMLEquations
  TYPE(OC_ComputationEnvironmentType) :: computationEnvironment
  TYPE(OC_ContextType) :: context
  TYPE(OC_ControlLoopType) :: controlLoop
  TYPE(OC_CoordinateSystemType) :: coordinateSystem
  TYPE(OC_DecompositionType) :: decomposition
  TYPE(OC_DecomposerType) :: decomposer
  TYPE(OC_EquationsType) :: equations
  TYPE(OC_EquationsSetType) :: equationsSet
  TYPE(OC_FieldType) :: geometricField,equationsSetField,dependentField,materialsField,sourceField
  TYPE(OC_FieldType) :: cellMLModelsField,cellMLStateField,cellMLIntermediateField,cellMLParametersField
  TYPE(OC_FieldsType) :: fields
  TYPE(OC_GeneratedMeshType) :: generatedMesh  
  TYPE(OC_MeshType) :: mesh
  TYPE(OC_ProblemType) :: problem
  TYPE(OC_RegionType) :: region,worldRegion
  TYPE(OC_SolverType) :: solver, linearSolver
  TYPE(OC_SolverEquationsType) :: solverEquations
  TYPE(OC_WorkGroupType) :: worldWorkGroup

  LOGICAL :: exportField
  
  !Generic CMISS variables
  
  INTEGER(OC_Intg) :: numberOfComputationalNodes,computationalNodeNumber
  INTEGER(OC_Intg) :: decompositionIndex,equationsSetIndex,cellMLIndex
  INTEGER(OC_Intg) :: err

  !-----------------------------------------------------------------------------------------------------------
  ! PROBLEM CONTROL PANEL
  !-----------------------------------------------------------------------------------------------------------

  !Intialise OpenCMISS
  CALL OC_Initialise(err)
  CALL OC_ErrorHandlingModeSet(OC_ERRORS_TRAP_ERROR,err)
  !Create a context
  CALL OC_Context_Initialise(context,err)
  CALL OC_Context_Create(CONTEXT_USER_NUMBER,context,err)
  CALL OC_ErrorHandlingModeSet(OC_ERRORS_TRAP_ERROR,err)
  CALL OC_Region_Initialise(worldRegion,err)
  CALL OC_Context_WorldRegionGet(context,worldRegion,err)
  
  !Get the computational nodes information
  CALL OC_ComputationEnvironment_Initialise(computationEnvironment,err)
  CALL OC_Context_ComputationEnvironmentGet(context,computationEnvironment,err)
  
  CALL OC_WorkGroup_Initialise(worldWorkGroup,err)
  CALL OC_ComputationEnvironment_WorldWorkGroupGet(computationEnvironment,worldWorkGroup,err)
  CALL OC_WorkGroup_NumberOfGroupNodesGet(worldWorkGroup,numberOfComputationalNodes,err)
  CALL OC_WorkGroup_GroupNodeNumberGet(worldWorkGroup,computationalNodeNumber,err)

  numberOfGlobalXElements=10

  !-----------------------------------------------------------------------------------------------------------
  ! COORDINATE SYSTEM
  !-----------------------------------------------------------------------------------------------------------  

  !Start the creation of a new RC coordinate system
  CALL OC_CoordinateSystem_Initialise(coordinateSystem,err)
  CALL OC_CoordinateSystem_CreateStart(COORDINATE_SYSTEM_USER_NUMBER,context,coordinateSystem,err)
  !Set the coordinate system to be 1D
  CALL OC_CoordinateSystem_DimensionSet(coordinateSystem,1,err)
  !Finish the creation of the coordinate system
  CALL OC_CoordinateSystem_CreateFinish(coordinateSystem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! REGION
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of the region
  CALL OC_Region_Initialise(region,err)
  CALL OC_Region_CreateStart(REGION_USER_NUMBER,worldRegion,region,err)
  CALL OC_Region_LabelSet(region,"cellml_split_reaction_diffusion_equation",err)
  !Set the regions coordinate system to the 1D RC coordinate system that we have created
  CALL OC_Region_CoordinateSystemSet(region,coordinateSystem,err)
  !Finish the creation of the region
  CALL OC_Region_CreateFinish(region,err)

  !-----------------------------------------------------------------------------------------------------------
  ! BASIS
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of a basis (default is linear lagrange)
  CALL OC_Basis_Initialise(basis,err)
  CALL OC_Basis_CreateStart(BASIS_USER_NUMBER,context,basis,err)
  !Set the basis to be a linear Lagrange basis
  CALL OC_Basis_NumberOfXiSet(basis,1,err)
  !Finish the creation of the basis
  CALL OC_Basis_CreateFinish(basis,err)

  !-----------------------------------------------------------------------------------------------------------
  ! MESH
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of a generated mesh in the region
  CALL OC_GeneratedMesh_Initialise(generatedMesh,err)
  CALL OC_GeneratedMesh_CreateStart(GENERATED_MESH_USER_NUMBER,region,generatedMesh,err)
  !Set up a regular 1D mesh
  CALL OC_GeneratedMesh_TypeSet(generatedMesh,OC_GENERATED_MESH_REGULAR_MESH_TYPE,err)
  !Set the default basis
  CALL OC_GeneratedMesh_BasisSet(generatedMesh,basis,err)   
  !Define the mesh on the region
  CALL OC_GeneratedMesh_ExtentSet(generatedMesh,[LENGTH],err)
  CALL OC_GeneratedMesh_NumberOfElementsSet(generatedMesh,[numberOfGlobalXElements],err)
  !Finish the creation of a generated mesh in the region
  CALL OC_Mesh_Initialise(mesh,err)
  CALL OC_GeneratedMesh_CreateFinish(generatedMesh,MESH_USER_NUMBER,mesh,err)

  !-----------------------------------------------------------------------------------------------------------
  ! DECOMPOSITION
  !-----------------------------------------------------------------------------------------------------------

  !Create a decomposition
  CALL OC_Decomposition_Initialise(decomposition,err)
  CALL OC_Decomposition_CreateStart(DECOMPOSITION_USER_NUMBER,mesh,decomposition,err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL OC_Decomposition_TypeSet(decomposition,OC_DECOMPOSITION_CALCULATED_TYPE,err)
  !Finish the decomposition
  CALL OC_Decomposition_CreateFinish(decomposition,err)

  !-----------------------------------------------------------------------------------------------------------
  ! DECOMPOSER
  !-----------------------------------------------------------------------------------------------------------

  CALL OC_Decomposer_Initialise(decomposer,err)
  CALL OC_Decomposer_CreateStart(DECOMPOSER_USER_NUMBER,region,worldWorkGroup,decomposer,err)
  !Add in the decomposition
  CALL OC_Decomposer_DecompositionAdd(decomposer,decomposition,decompositionIndex,err)
  !Finish the decomposer
  CALL OC_Decomposer_CreateFinish(decomposer,err)
  
  !-----------------------------------------------------------------------------------------------------------
  ! GEOMETRIC FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Start to create a default (geometric) field on the region
  CALL OC_Field_Initialise(geometricField,err)
  CALL OC_Field_CreateStart(GEOMETRIC_FIELD_USER_NUMBER,region,geometricField,err)
  !Set the decomposition to use
  CALL OC_Field_DecompositionSet(geometricField,decomposition,err)
  !Set the domain to be used by the field components.
  CALL OC_Field_ComponentMeshComponentSet(geometricField,OC_FIELD_U_VARIABLE_TYPE,1,1,err)
  !Finish creating the field
  CALL OC_Field_CreateFinish(geometricField,err)

  !Update the geometric field parameters
  CALL OC_GeneratedMesh_GeometricParametersCalculate(generatedMesh,geometricField,err)

  !-----------------------------------------------------------------------------------------------------------
  ! EQUATIONS SETS
  !-----------------------------------------------------------------------------------------------------------

  !Create the cellml reaction with split reaction diffusion equations_set 
  CALL OC_EquationsSet_Initialise(equationsSet,err)
  CALL OC_Field_Initialise(equationsSetField,err)
  CALL OC_EquationsSet_CreateStart(EQUATIONS_SET_USER_NUMBER,region,geometricField,[OC_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & OC_EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE,OC_EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE], &
    & EQUATIONS_SET_FIELD_USER_NUMBER,equationsSetField,equationsSet,err)
  !Finish creating the equations set
  CALL OC_EquationsSet_CreateFinish(equationsSet,err)

  !-----------------------------------------------------------------------------------------------------------
  ! DEPENDENT FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set dependent field variables (primary and secondary variable)
  CALL OC_Field_Initialise(dependentField,err)
  CALL OC_EquationsSet_DependentCreateStart(equationsSet,DEPENDENT_FIELD_USER_NUMBER,dependentField,err)
  CALL OC_Field_VariableLabelSet(dependentField,OC_FIELD_U_VARIABLE_TYPE,"U",err)
  CALL OC_Field_VariableLabelSet(dependentField,OC_FIELD_DELUDELN_VARIABLE_TYPE,"DELUDELN",err)
  !Finish the equations set dependent field variables
  CALL OC_EquationsSet_DependentCreateFinish(equationsSet,err)

  !-----------------------------------------------------------------------------------------------------------
  ! MATERIAL FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set material field variables
  !by default 2 components for reaction-diffusion i.e. diffusion coefficient in x direction 
  !set constant spatially = 1, and storage coefficient set to 1
  CALL OC_Field_Initialise(materialsField,err)
  CALL OC_EquationsSet_MaterialsCreateStart(equationsSet,MATERIALS_FIELD_USER_NUMBER,materialsField,err)
  CALL OC_Field_VariableLabelSet(materialsField,OC_FIELD_U_VARIABLE_TYPE,"Material",err)
  !Finish the equations set materials field variables
  CALL OC_EquationsSet_MaterialsCreateFinish(equationsSet,err)
  CALL OC_Field_ComponentValuesInitialise(materialsField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
   & 1,0.5_OC_RP,err) !diffusion coefficent in x
  CALL OC_Field_ComponentValuesInitialise(materialsField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
   & 2,1.0_OC_RP,err) ! storage coefficient
  
  !-----------------------------------------------------------------------------------------------------------
  ! SOURCE FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Set up source field for reaction diffusion equation set. Note that for the split problem subtype, 
  !the source field is not used at all.
  CALL OC_Field_Initialise(sourceField,err)
  CALL OC_EquationsSet_SourceCreateStart(equationsSet,SOURCE_FIELD_USER_NUMBER,sourceField,err)
  CALL OC_Field_VariableLabelSet(sourceField,OC_FIELD_U_VARIABLE_TYPE,"Source",err)
  !Finish the equations set source field variables
  CALL OC_EquationsSet_SourceCreateFinish(equationsSet,err)
  CALL OC_Field_ComponentValuesInitialise(sourceField,OC_FIELD_U_VARIABLE_TYPE, &
    & OC_FIELD_VALUES_SET_TYPE,1,0.0_OC_RP,err)

  !-----------------------------------------------------------------------------------------------------------
  ! CELLML FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Start to set up CellML Fields

  !Create the CellML environment
  CALL OC_CellML_Initialise(cellML,err)
  CALL OC_CellML_CreateStart(CELLML_USER_NUMBER,region,cellML,err)
  !Import a constant source (i.e. rate of generation/depletion is zero) model from a file
  CALL OC_CellML_ModelImport(cellML,"constant_rate.xml",constantModelIndex,err)
  !Speify the variables in the imported model that will be used. 
  
  ! Now we have imported all the models we are able to specify which variables from the model we want:
  ! - to set from this side
  !These are effectively parameters that won't change in the course of the ode solving for one time step. 
  !i.e. fixed before running cellml, known in opencmiss and changed only in opencmiss - components of the parameters field
  CALL OC_CellML_VariableSetAsKnown(cellML,constantModelIndex,"dude/param",err)
  ! - to get from the CellML side. variables in cellml model that are not state variables, but are dependent on 
  !independent and state variables. - components of intermediate field
  CALL OC_CellML_VariableSetAsWanted(cellML,constantModelIndex,"dude/intmd",err)
  !Finish the CellML environment
  CALL OC_CellML_CreateFinish(cellML,err)

  !Start the creation of CellML <--> OpenCMISS field maps
  CALL OC_CellML_FieldMapsCreateStart(cellML,err)
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
 
  CALL OC_CellML_CreateFieldToCellMLMap(cellML,dependentField,OC_FIELD_U_VARIABLE_TYPE,1,OC_FIELD_VALUES_SET_TYPE, &
    & constantModelIndex,"dude/ca",OC_FIELD_VALUES_SET_TYPE,err)
  CALL OC_CellML_CreateCellMLToFieldMap(cellML,constantModelIndex,"dude/ca",OC_FIELD_VALUES_SET_TYPE, &
    & dependentField,OC_FIELD_U_VARIABLE_TYPE,1,OC_FIELD_VALUES_SET_TYPE,err)
  !Finish the creation of CellML <--> OpenCMISS field maps
  CALL OC_CellML_FieldMapsCreateFinish(cellML,err)

  !set initial value of the dependent field/state variable, Ca concentration.
  CALL OC_Field_ComponentValuesInitialise(dependentField,OC_FIELD_U_VARIABLE_TYPE, &
    & OC_FIELD_VALUES_SET_TYPE,1,0.0_OC_RP,err)
  nodeIdx=2
  CALL OC_Decomposition_NodeDomainGet(decomposition,nodeIdx,1,nodeDomain,err)
  IF(nodeDomain==computationalNodeNumber) THEN
    CALL OC_Field_ParameterSetUpdateNode(dependentField,OC_FIELD_U_VARIABLE_TYPE, &
     & OC_FIELD_VALUES_SET_TYPE, &
     & 1,1,nodeIdx,1,0.0_OC_RP,err) 
  ENDIF
  !Start the creation of the CellML models field. This field is an integer field that stores which nodes have which cellml 
  !model
  CALL OC_Field_Initialise(cellMLModelsField,err)
  CALL OC_CellML_ModelsFieldCreateStart(cellML, CELLML_MODELS_FIELD_USER_NUMBER, &
    & cellMLModelsField,err)
  !Finish the creation of the CellML models field
  CALL OC_CellML_ModelsFieldCreateFinish(cellML,err)
  !The CellMLModelsField is an integer field that stores which model is being used by which node.
  !By default all field parameters have default model value of 1, i.e. the first model. But, this command below is 
  !for example purposes
  CALL OC_Field_ComponentValuesInitialise(cellMLModelsField,OC_FIELD_U_VARIABLE_TYPE, &
    & OC_FIELD_VALUES_SET_TYPE,1,1_OC_Intg,err)

  !Set up the models field
  !DO N=1,(numberOfGlobalXElements+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)*(NUMBER_GLOBAL_Z_ELEMENTS+1)
  !  IF(N < 5) THEN
  !    CELL_TYPE = 1
  !  ELSE
  !    CELL_TYPE = 2
  !  ENDIF
  !  CALL OC_FieldParameterSetUpdateNode(cellMLModelsField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,N,1,CELL_TYPE,err)
  !END DO
  !CALL OC_FieldParameterSetUpdateStart(cellMLModelsField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,err)
  !CALL OC_FieldParameterSetUpdateFinish(cellMLModelsField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,err)

  !Start the creation of the cellML state field
  CALL OC_Field_Initialise(cellMLStateField,err)
  CALL OC_CellML_StateFieldCreateStart(cellML,CELLML_STATE_FIELD_USER_NUMBER,cellMLStateField,err)
  !Finish the creation of the CellML state field
  CALL OC_CellML_StateFieldCreateFinish(cellML,err)

  !Start the creation of the CellML intermediate field
  CALL OC_Field_Initialise(cellMLIntermediateField,err)
  CALL OC_CellML_IntermediateFieldCreateStart(cellML,CELLML_INTERMEDIATE_FIELD_USER_NUMBER,cellMLIntermediateField,err)
  !Finish the creation of the CellML intermediate field
  CALL OC_CellML_IntermediateFieldCreateFinish(cellML,err)

  !Start the creation of CellML parameters field
  CALL OC_Field_Initialise(cellMLParametersField,err)
  CALL OC_CellML_ParametersFieldCreateStart(cellML,CELLML_PARAMETERS_FIELD_USER_NUMBER,cellMLParametersField,err)
  !Finish the creation of CellML parameters
  CALL OC_CellML_ParametersFieldCreateFinish(cellML,err)

  !-----------------------------------------------------------------------------------------------------------
  ! EQUATIONS SETS
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set equations
  CALL OC_Equations_Initialise(equations,err)
  CALL OC_EquationsSet_EquationsCreateStart(equationsSet,equations,err)
  !Set the equations matrices sparsity type
  CALL OC_Equations_SparsityTypeSet(equations,OC_EQUATIONS_SPARSE_MATRICES,err)
  !Set the equations set output
  CALL OC_Equations_OutputTypeSet(equations,OC_EQUATIONS_NO_OUTPUT,err)
  !CALL OC_EquationsOutputTypeSet(equations,OC_EQUATIONS_TIMING_OUTPUT,err)
  !CALL OC_EquationsOutputTypeSet(equations,OC_EQUATIONS_MATRIX_OUTPUT,err)
  !CALL OC_equationsOutputTypeSet(equations,OC_EQUATIONS_ELEMENT_MATRIX_OUTPUT,err)
  !Finish the equations set equations
  CALL OC_EquationsSet_EquationsCreateFinish(equationsSet,err)

  !-----------------------------------------------------------------------------------------------------------
  ! PROBLEM
  !-----------------------------------------------------------------------------------------------------------

  !Create the problem
  CALL OC_Problem_Initialise(problem,err)
  CALL OC_Problem_CreateStart(PROBLEM_USER_NUMBER,context,[OC_PROBLEM_CLASSICAL_FIELD_CLASS, &
    & OC_PROBLEM_REACTION_DIFFUSION_EQUATION_TYPE,OC_PROBLEM_CELLML_REAC_INTEG_REAC_DIFF_STRANG_SPLIT_SUBTYPE],problem,err)
  !Finish the creation of a problem.
  CALL OC_Problem_CreateFinish(problem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! PROBLEM CONTROL
  !-----------------------------------------------------------------------------------------------------------

  !Create the problem control
  CALL OC_Problem_ControlLoopCreateStart(problem,err)
  !Get the control loop
  CALL OC_ControlLoop_Initialise(controlLoop,err)
  CALL OC_Problem_ControlLoopGet(problem,OC_CONTROL_LOOP_NODE,controlLoop,err)
  !Set the times
  CALL OC_ControlLoop_TimesSet(controlLoop,0.0_OC_RP,0.5_OC_RP,0.01_OC_RP,err)
  CALL OC_ControlLoop_OutputTypeSet(controlLoop,OC_CONTROL_LOOP_PROGRESS_OUTPUT,err)
  !Finish creating the problem control loop
  CALL OC_Problem_ControlLoopCreateFinish(problem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! PROBLEM SOLVERS
  !-----------------------------------------------------------------------------------------------------------

  !Set up the problem solvers for Strang splitting
  CALL OC_Problem_SolversCreateStart(problem,err)
  !First solver is a DAE solver
  CALL OC_Solver_Initialise(solver,err)
  CALL OC_Problem_SolverGet(problem,OC_CONTROL_LOOP_NODE,1,solver,err)
  CALL OC_Solver_DAESolverTypeSet(solver,OC_SOLVER_DAE_EULER,err)
  CALL OC_Solver_DAETimeStepSet(solver,0.0000001_OC_RP,err)
  CALL OC_Solver_OutputTypeSet(solver,OC_SOLVER_MATRIX_OUTPUT,err)

  !Second solver is the dynamic solver for solving the parabolic equation
  CALL OC_Solver_Initialise(solver,err)
  CALL OC_Solver_Initialise(linearSolver,err)
  CALL OC_Problem_SolverGet(problem,OC_CONTROL_LOOP_NODE,2,solver,err)
  !set theta - backward vs forward time step parameter
  CALL OC_Solver_DynamicThetaSet(solver,1.0_OC_RP,err)
  !CALL OC_SolverOutputTypeSet(solver,OC_SOLVER_NO_OUTPUT,err)
  CALL OC_Solver_OutputTypeSet(solver,OC_SOLVER_TIMING_OUTPUT,err)
  !get the dynamic linear solver from the solver
  CALL OC_Solver_DynamicLinearSolverGet(solver,linearSolver,err)
  CALL OC_Solver_LibraryTypeSet(linearSolver,OC_SOLVER_LAPACK_LIBRARY,err)
  CALL OC_Solver_LinearTypeSet(linearSolver,OC_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,err)
  CALL OC_Solver_LinearDirectTypeSet(linearSolver,OC_SOLVER_DIRECT_LU,err)
  !CALL OC_SolverLibraryTypeSet(linearSolver,OC_SOLVER_OC__LIBRARY,err)
  !CALL OC_SolverLinearTypeSet(linearSolver,OC_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,err)
  !CALL OC_SolverLibraryTypeSet(linearSolver,OC_SOLVER_MUMPS_LIBRARY,err)
  !CALL OC_Solver_LinearIterativeMaximumIterationsSet(linearSolver,10000,err)

  !Third solver is another DAE solver
  CALL OC_Solver_Initialise(solver,err)
  CALL OC_Problem_SolverGet(problem,OC_CONTROL_LOOP_NODE,3,solver,err)
  CALL OC_Solver_DAESolverTypeSet(solver,OC_SOLVER_DAE_EULER,err)
  CALL OC_Solver_DAETimeStepSet(solver,0.0000001_OC_RP,err)
  CALL OC_Solver_OutputTypeSet(solver,OC_SOLVER_MATRIX_OUTPUT,err)
  !CALL OC_SolverOutputTypeSet(solver,OC_SOLVER_TIMING_OUTPUT,err)
  !CALL OC_SolverOutputTypeSet(solver,OC_SOLVER_SOLVER_OUTPUT,err)
  !CALL OC_SolverOutputTypeSet(solver,OC_SOLVER_PROGRESS_OUTPUT,err)
  !Finish the creation of the problem solver
  CALL OC_Problem_SolversCreateFinish(problem,err)

  !Start the creation of the problem solver CellML equations
  CALL OC_Problem_CellMLEquationsCreateStart(problem,err)
  !Get the first solver  
  !Get the CellML equations
  CALL OC_Solver_Initialise(solver,err)
  CALL OC_Problem_SolverGet(problem,OC_CONTROL_LOOP_NODE,1,solver,err)
  CALL OC_CellMLEquations_Initialise(cellMLEquations,err)
  CALL OC_Solver_CellMLEquationsGet(solver,cellMLEquations,err)
  !Add in the CellML environement
  CALL OC_CellMLEquations_CellMLAdd(cellMLEquations,cellML,cellMLIndex,err)

  !Get the third solver  
  !Get the CellML equations
  CALL OC_Solver_Initialise(solver,err)
  CALL OC_Problem_SolverGet(problem,OC_CONTROL_LOOP_NODE,3,solver,err)
  CALL OC_CellMLEquations_Initialise(cellMLEquations,err)
  CALL OC_Solver_CellMLEquationsGet(solver,cellMLEquations,err)
  !Add in the CellML environement
  CALL OC_CellMLEquations_CellMLAdd(cellMLEquations,cellML,cellMLIndex,err)

  !Finish the creation of the problem solver CellML equations
  CALL OC_Problem_CellMLEquationsCreateFinish(problem,err)

  !Start the creation of the problem solver equations
  CALL OC_Problem_SolverEquationsCreateStart(problem,err)
  !Get the second solver  
  !Get the solver equations
  CALL OC_Solver_Initialise(solver,err)
  CALL OC_Problem_SolverGet(problem,OC_CONTROL_LOOP_NODE,2,solver,err)
  CALL OC_SolverEquations_Initialise(solverEquations,err)
  CALL OC_Solver_SolverEquationsGet(solver,solverEquations,err)
  !Set the solver equations sparsity
  !CALL OC_SolverEquationsSparsityTypeSet(solverEquations,OC_SOLVER_EQUATIONS_SPARSE_MATRICES,err)
  CALL OC_SolverEquations_SparsityTypeSet(solverEquations,OC_SOLVER_SPARSE_MATRICES,err)  
  !Add in the equations set
  CALL OC_SolverEquations_EquationsSetAdd(solverEquations,equationsSet,equationsSetIndex,err)
  !Finish the creation of the problem solver equations
  CALL OC_Problem_SolverEquationsCreateFinish(problem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! BOUNDARY CONDITIONS
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set boundary conditions
  boundaryConditionNodes = [1,11]
  CALL OC_BoundaryConditions_Initialise(boundaryConditions,err)
  CALL OC_SolverEquations_BoundaryConditionsCreateStart(solverEquations,boundaryConditions,err)
  DO nodeIdx=1,2
    nodeNumber = boundaryConditionNodes(nodeIdx)
    CALL OC_Decomposition_NodeDomainGet(decomposition,nodeNumber,1,nodeDomain,err)
    IF(nodeDomain==computationalNodeNumber) THEN
      CALL OC_BoundaryConditions_SetNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE,1,OC_NO_GLOBAL_DERIV, &
        & nodeNumber,1,OC_BOUNDARY_CONDITION_FIXED,1.5_OC_RP,err)
      !Need to set CellML model to zero at the nodes at which value of ca has been fixed.
      CALL OC_Field_ParameterSetUpdateNode(cellMLModelsField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,&
        & 1,1,nodeNumber,1,0_OC_Intg,err) 
    ENDIF
  ENDDO
  CALL OC_SolverEquations_BoundaryConditionsCreateFinish(solverEquations,err)

  !-----------------------------------------------------------------------------------------------------------
  ! SOLVER
  !-----------------------------------------------------------------------------------------------------------

  !Solve the problem
  CALL OC_Problem_Solve(problem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! OUTPUT
  !-----------------------------------------------------------------------------------------------------------
 
  exportField=.TRUE.
  IF(exportField) THEN
    CALL OC_Fields_Initialise(fields,err)
    CALL OC_Fields_Create(region,fields,err)
    CALL OC_Fields_NodesExport(fields,"cellml_split_reaction_diffusion_equation","FORTRAN",err)
    CALL OC_Fields_ElementsExport(fields,"cellml_split_reaction_diffusion_equation","FORTRAN",err)
    CALL OC_Fields_Finalise(fields,err)
  ENDIF

  !Destroy the context
  CALL OC_Context_Destroy(context,err)
  !Finalise OpenCMISS
  CALL OC_Finalise(err)
  
  WRITE(*,'(A)') "Program successfully completed."

  STOP
  
END PROGRAM CellMLSplitReactionDiffusionEquation
