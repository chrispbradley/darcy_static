PROGRAM DarcyStaticExample

  USE OpenCMISS
 
  IMPLICIT NONE

  !-----------------------------------------------------------------------------------------------------------
  ! PROGRAM VARIABLES AND TYPES
  !-----------------------------------------------------------------------------------------------------------

  !Test program parameters
  REAL(OC_RP), PARAMETER :: HEIGHT=1.0_OC_RP
  REAL(OC_RP), PARAMETER :: WIDTH=1.0_OC_RP
  REAL(OC_RP), PARAMETER :: LENGTH=1.0_OC_RP

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
  INTEGER(OC_Intg), PARAMETER :: ANALYTIC_FIELD_USER_NUMBER=13
  INTEGER(OC_Intg), PARAMETER :: EQUATIONS_SET_USER_NUMBER=14
  INTEGER(OC_Intg), PARAMETER :: PROBLEM_USER_NUMBER=15

  !Program variables
  INTEGER(OC_Intg) :: numberOfGlobalXElements,numberOfGlobalYElements,numberOfGlobalZElements
  INTEGER(OC_Intg) :: numberOfDimensions
  INTEGER(OC_Intg) :: basisType
  INTEGER(OC_Intg) :: numberOfGaussXi
  INTEGER(OC_Intg) :: basisXiInterpolation
  INTEGER(OC_Intg) :: maximumIterations
  INTEGER(OC_Intg) :: restartValue
  INTEGER(OC_Intg) :: numberOfFixedWallNodes
  INTEGER(OC_Intg) :: numberOfInletWallNodes
  INTEGER(OC_Intg) :: equationsDarcyOutput
  INTEGER(OC_Intg) :: componentIdx
  INTEGER(OC_Intg) :: nodeNumber
  INTEGER(OC_Intg) :: nodeIdx
  INTEGER(OC_Intg) :: condition
  INTEGER(OC_Intg) :: linearSolverDarcyOutputType
  INTEGER(OC_Intg) :: analyticalType
  INTEGER, ALLOCATABLE :: fixedWallNodes(:)
  INTEGER, ALLOCATABLE :: inletWallNodes(:)
  REAL(OC_RP) :: initialFieldValue(3)
  REAL(OC_RP) :: boundaryConditionValues(3)
  REAL(OC_RP) :: divergenceTolerance
  REAL(OC_RP) :: relativeTolerance
  REAL(OC_RP) :: absoluteTolerance
  REAL(OC_RP) :: linesearchAlpha
  REAL(OC_RP) :: boundaryConditionValue
  REAL(OC_RP) :: porosityParameter,permOverVisParameter
  LOGICAL :: directoryExists
  LOGICAL :: exportFieldIO
  LOGICAL :: linearSolverDarcyDirectFlag
  LOGICAL :: fixedWallNodesFlag
  LOGICAL :: inletWallNodesFlag

  !CMISS variables
  TYPE(OC_BasisType) :: basis
  TYPE(OC_BoundaryConditionsType) :: boundaryConditionsDarcy
  TYPE(OC_ComputationEnvironmentType) :: computationEnvironment
  TYPE(OC_ContextType) :: context
  TYPE(OC_ControlLoopType) :: controlLoop
  TYPE(OC_CoordinateSystemType) :: coordinateSystem
  TYPE(OC_DecomposerType) :: decomposer
  TYPE(OC_DecompositionType) :: decomposition
  TYPE(OC_EquationsType) :: equationsDarcy
  TYPE(OC_EquationsSetType) :: equationsSetDarcy
  TYPE(OC_FieldType) :: geometricField,analyticField
  TYPE(OC_FieldType) :: equationsSetField
  TYPE(OC_FieldType) :: dependentFieldDarcy
  TYPE(OC_FieldType) :: materialsFieldDarcy
  TYPE(OC_FieldsType) :: fields
  TYPE(OC_GeneratedMeshType) :: generatedMesh
  TYPE(OC_MeshType) :: mesh
  TYPE(OC_ProblemType) :: problem
  TYPE(OC_SolverType) :: linearSolverDarcy
  TYPE(OC_SolverEquationsType) :: solverEquationsDarcy
  TYPE(OC_RegionType) :: region
  TYPE(OC_RegionType) :: WorldRegion
  TYPE(OC_WorkGroupType) :: worldWorkGroup

  !Generic CMISS variables
  INTEGER(OC_Intg) :: numberOfComputationalNodes,computationalNodeNumber
  INTEGER(OC_Intg) :: decompositionIndex,equationsSetIndex,boundaryNodeDomain,err

  !Intialise OpenCMISS
  CALL OC_Initialise(err)
  CALL OC_ErrorHandlingModeSet(OC_ERRORS_TRAP_ERROR,err)
  !Create a context
  CALL OC_Context_Initialise(context,err)
  CALL OC_Context_Create(CONTEXT_USER_NUMBER,context,err)
  CALL OC_Region_Initialise(worldRegion,err)
  CALL OC_Context_WorldRegionGet(context,worldRegion,err)
 
  !Get the computational nodes information
  CALL OC_ComputationEnvironment_Initialise(computationEnvironment,err)
  CALL OC_Context_ComputationEnvironmentGet(context,computationEnvironment,err)
  
  CALL OC_WorkGroup_Initialise(worldWorkGroup,err)
  CALL OC_ComputationEnvironment_WorldWorkGroupGet(computationEnvironment,worldWorkGroup,err)
  CALL OC_WorkGroup_NumberOfGroupNodesGet(worldWorkGroup,numberOfComputationalNodes,err)
  CALL OC_WorkGroup_GroupNodeNumberGet(worldWorkGroup,computationalNodeNumber,err)

  !-----------------------------------------------------------------------------------------------------------
  ! PROBLEM CONTROL PANEL
  !-----------------------------------------------------------------------------------------------------------

  numberOfGlobalXElements=3
  numberOfGlobalYElements=3
  numberOfGlobalZElements=3
  numberOfDimensions=3
  basisType=OC_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  basisXiInterpolation=OC_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  !Set initial values
  initialFieldValue(1)=0.0_OC_RP
  initialFieldValue(2)=0.0_OC_RP
  initialFieldValue(3)=0.0_OC_RP
  !Set boundary conditions
  fixedWallNodesFlag=.TRUE.
  inletWallNodesFlag=.TRUE.
  IF(fixedWallNodesFlag) THEN
    numberOfFixedWallNodes=48
    ALLOCATE(fixedWallNodes(numberOfFixedWallNodes))
    fixedWallNodes=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,17,20,21,24,25, &
      & 28,29,32,33,36,37,40,41,44,45,48]
  ENDIF
  IF(inletWallNodesFlag) THEN
    numberOfInletWallNodes=4
    ALLOCATE(inletWallNodes(numberOfInletWallNodes))
    inletWallNodes=[18,19,34,35]
    !Set initial boundary conditions
    boundaryConditionValues(1)=0.0_OC_RP
    boundaryConditionValues(2)=1.0_OC_RP
    boundaryConditionValues(3)=0.0_OC_RP
  ENDIF
  !Set material parameters
  porosityParameter=0.5_OC_RP
  permOverVisParameter=0.8_OC_RP
  !Set number of Gauss points
  numberOfGaussXi=3
  !Set output parameter
  !(NoOutput/ProgressOutput/TimingOutput/SolverOutput/SolverMatrixOutput)
  linearSolverDarcyOutputType=OC_SOLVER_NO_OUTPUT
  !(NoOutput/TimingOutput/MatrixOutput/ElementOutput)
  equationsDarcyOutput=OC_EQUATIONS_NO_OUTPUT
  !Set solver parameters
  linearSolverDarcyDirectFlag=.FALSE.
  relativeTolerance=1.0E-10_OC_RP !default: 1.0E-05_OC_RP
  absoluteTolerance=1.0E-10_OC_RP !default: 1.0E-10_OC_RP
  divergenceTolerance=1.0E5_OC_RP !default: 1.0E5
  maximumIterations=10000_OC_Intg !default: 100000
  restartValue=3000_OC_Intg !default: 30
  linesearchAlpha=1.0_OC_RP

  !IF(numberOfDimensions==2) analyticalType=OC_EQUATIONS_SET_DARCY_EQUATION_TWO_DIM_2
  IF(numberOfDimensions==2) analyticalType=OC_EQUATIONS_SET_DARCY_EQUATION_TWO_DIM_3
  !IF(numberOfDimensions==3) analyticalType=OC_EQUATIONS_SET_DARCY_EQUATION_THREE_DIM_2
  IF(numberOfDimensions==3) analyticalType=OC_EQUATIONS_SET_DARCY_EQUATION_THREE_DIM_3

  !-----------------------------------------------------------------------------------------------------------
  ! COORDINATE SYSTEM
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of a new RC coordinate system
  CALL OC_CoordinateSystem_Initialise(coordinateSystem,err)
  CALL OC_CoordinateSystem_CreateStart(COORDINATE_SYSTEM_USER_NUMBER,context,coordinateSystem,err)
  !Set the coordinate system dimension
  CALL OC_CoordinateSystem_DimensionSet(coordinateSystem,numberOfDimensions,err)
  !Finish the creation of the coordinate system
  CALL OC_CoordinateSystem_CreateFinish(coordinateSystem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! REGION
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of a new region
  CALL OC_Region_Initialise(region,err)
  CALL OC_Region_CreateStart(REGION_USER_NUMBER,worldRegion,region,err)
  CALL OC_Region_LabelSet(region,"DarcyRegion",err)
  !Set the regions coordinate system as defined above
  CALL OC_Region_CoordinateSystemSet(region,coordinateSystem,err)
  !Finish the creation of the region
  CALL OC_Region_CreateFinish(region,err)

  !-----------------------------------------------------------------------------------------------------------
  ! BASIS
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of new bases
  CALL OC_Basis_Initialise(basis,err)
  CALL OC_Basis_CreateStart(BASIS_USER_NUMBER,context,basis,err)
  !Set the basis type (Lagrange/Simplex)
  CALL OC_Basis_TypeSet(basis,basisType,err)
  !Set the basis xi number
  CALL OC_Basis_NumberOfXiSet(basis,numberOfDimensions,err)
  !Set the basis xi interpolation and number of Gauss points
  IF(numberOfDimensions==2) THEN
    CALL OC_Basis_InterpolationXiSet(basis,[basisXiInterpolation,basisXiInterpolation],err)
    CALL OC_Basis_QuadratureNumberOfGaussXiSet(basis,[numberOfGaussXi,numberOfGaussXi],err)
  ELSE
    CALL OC_Basis_InterpolationXiSet(basis,[basisXiInterpolation,basisXiInterpolation,basisXiInterpolation],err)
    CALL OC_Basis_QuadratureNumberOfGaussXiSet(basis,[numberOfGaussXi,numberOfGaussXi,numberOfGaussXi],err)
    CALL OC_Basis_QuadratureLocalFaceGaussEvaluateSet(basis,.TRUE.,err)
  ENDIF
  !Finish the creation of the basis
  CALL OC_Basis_CreateFinish(basis,err)

  !-----------------------------------------------------------------------------------------------------------
  ! MESH
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of a generated mesh in the region
  CALL OC_GeneratedMesh_Initialise(generatedMesh,err)
  CALL OC_GeneratedMesh_CreateStart(GENERATED_MESH_USER_NUMBER,region,generatedMesh,err)
  !Set up a regular x*y*z mesh
  CALL OC_GeneratedMesh_TypeSet(generatedMesh,OC_GENERATED_MESH_REGULAR_MESH_TYPE,err)
  !Set the default basis
  CALL OC_GeneratedMesh_BasisSet(generatedMesh,[basis],err)
  !Define the mesh on the region
  IF(numberOfDimensions==2) THEN
     CALL OC_GeneratedMesh_ExtentSet(generatedMesh,[WIDTH,HEIGHT],err)
     CALL OC_GeneratedMesh_NumberOfElementsSet(generatedMesh,[numberOfGlobalXElements,numberOfGlobalYElements],err)
  ELSE
     CALL OC_GeneratedMesh_ExtentSet(generatedMesh,[WIDTH,HEIGHT,LENGTH],err)
     CALL OC_GeneratedMesh_NumberOfElementsSet(generatedMesh,[numberOfGlobalXElements,numberOfGlobalYElements, &
       & numberOfGlobalZElements],err)
  ENDIF
  !Finish the creation of a generated mesh in the region
  CALL OC_Mesh_Initialise(mesh,err)
  CALL OC_GeneratedMesh_CreateFinish(generatedMesh,MESH_USER_NUMBER,mesh,err)

  !-----------------------------------------------------------------------------------------------------------
  ! MESH DECOMPOSITION
  !-----------------------------------------------------------------------------------------------------------

  !Create a decomposition
  CALL OC_Decomposition_Initialise(decomposition,err)
  CALL OC_Decomposition_CreateStart(DECOMPOSITION_USER_NUMBER,mesh,decomposition,err)
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
  !Set the field type
  CALL OC_Field_TypeSet(geometricField,OC_FIELD_GEOMETRIC_TYPE,err)
  !Set the decomposition to use
  CALL OC_Field_DecompositionSet(geometricField,decomposition,err)
  !Set the scaling to use
  CALL OC_Field_ScalingTypeSet(geometricField,OC_FIELD_NO_SCALING,err)
  !Set the mesh component to be used by the field components.
  DO componentIdx=1,numberOfDimensions
     CALL OC_Field_ComponentMeshComponentSet(geometricField,OC_FIELD_U_VARIABLE_TYPE,componentIdx, &
       & 1,err)
  ENDDO
  !Finish creating the field
  CALL OC_Field_CreateFinish(geometricField,err)
  !Update the geometric field parameters
  CALL OC_GeneratedMesh_GeometricParametersCalculate(generatedMesh,geometricField,err)

  !-----------------------------------------------------------------------------------------------------------
  ! EQUATIONS SETS
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set for Static Darcy
  CALL OC_EquationsSet_Initialise(equationsSetDarcy,err)
  CALL OC_Field_Initialise(equationsSetField,err)
  CALL OC_EquationsSet_CreateStart(EQUATIONS_SET_USER_NUMBER,region,geometricField,[OC_EQUATIONS_SET_FLUID_MECHANICS_CLASS, &
    & OC_EQUATIONS_SET_DARCY_EQUATION_TYPE,OC_EQUATIONS_SET_STANDARD_DARCY_SUBTYPE],EQUATIONS_SET_FIELD_USER_NUMBER, &
    & equationsSetField,equationsSetDarcy,err)
  !Finish creating the equations set
  CALL OC_EquationsSet_CreateFinish(equationsSetDarcy,err)

  !-----------------------------------------------------------------------------------------------------------
  ! DEPENDENT FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set dependent field variables for Static Darcy
  CALL OC_Field_Initialise(dependentFieldDarcy,err)
  CALL OC_EquationsSet_DependentCreateStart(equationsSetDarcy,DEPENDENT_FIELD_USER_NUMBER,dependentFieldDarcy,err)
  !Set the mesh component to be used by the field components.
  DO componentIdx=1,numberOfDimensions
    CALL OC_Field_ComponentMeshComponentSet(dependentFieldDarcy,OC_FIELD_U_VARIABLE_TYPE,componentIdx,1,err)
    CALL OC_Field_ComponentMeshComponentSet(dependentFieldDarcy,OC_FIELD_DELUDELN_VARIABLE_TYPE,componentIdx,1,err)
  ENDDO
  componentIdx=numberOfDimensions+1
  CALL OC_Field_ComponentMeshComponentSet(dependentFieldDarcy,OC_FIELD_U_VARIABLE_TYPE,componentIdx,1,err)
  CALL OC_Field_ComponentMeshComponentSet(dependentFieldDarcy,OC_FIELD_DELUDELN_VARIABLE_TYPE,componentIdx,1,err)
  !Finish the equations set dependent field variables
  CALL OC_EquationsSet_DependentCreateFinish(equationsSetDarcy,err)

  !Initialise dependent field (velocity components)
  DO componentIdx=1,numberOfDimensions
    CALL OC_Field_ComponentValuesInitialise(dependentFieldDarcy,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
      & componentIdx,initialFieldValue(componentIdx),err)
  ENDDO
  
  !-----------------------------------------------------------------------------------------------------------
  ! MATERIALS FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set materials field variables for Static Darcy
  CALL OC_Field_Initialise(materialsFieldDarcy,err)
  CALL OC_EquationsSet_MaterialsCreateStart(equationsSetDarcy,MATERIALS_FIELD_USER_NUMBER,materialsFieldDarcy,err)
  !Finish the equations set materials field variables
  CALL OC_EquationsSet_MaterialsCreateFinish(equationsSetDarcy,err)
  CALL OC_Field_ComponentValuesInitialise(materialsFieldDarcy,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
    & 1,porosityParameter,err)
  CALL OC_Field_ComponentValuesInitialise(materialsFieldDarcy,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
    & 2,permOverVisParameter,err)

  !-----------------------------------------------------------------------------------------------------------
  ! ANALYTIC FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set analytic field variables for static Darcy
  !CALL OC_Field_Initialise(analyticField,err)
  !CALL OC_EquationsSet_AnalyticCreateStart(equationsSetDarcy,analyticalType,ANALYTIC_FIELD_USER_NUMBER,analyticField,err)
  !Finish the equations set analytic field variables
  !CALL OC_EquationsSet_AnalyticCreateFinish(equationsSetDarcy,err)

  !-----------------------------------------------------------------------------------------------------------
  ! EQUATIONS
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set equations
  CALL OC_Equations_Initialise(equationsDarcy,err)
  CALL OC_EquationsSet_EquationsCreateStart(equationsSetDarcy,equationsDarcy,err)
  !Set the equations matrices sparsity type
  CALL OC_Equations_SparsityTypeSet(equationsDarcy,OC_EQUATIONS_SPARSE_MATRICES,err)
  !Set the equations set output
  CALL OC_Equations_OutputTypeSet(equationsDarcy,equationsDarcyOutput,err)
  !Finish the equations set equations
  CALL OC_EquationsSet_EquationsCreateFinish(equationsSetDarcy,err)

  !-----------------------------------------------------------------------------------------------------------
  ! PROBLEM
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of a problem.
  CALL OC_Problem_Initialise(problem,err)
  CALL OC_ControlLoop_Initialise(controlLoop,err)
  CALL OC_Problem_CreateStart(PROBLEM_USER_NUMBER,context,[OC_PROBLEM_FLUID_MECHANICS_CLASS,OC_PROBLEM_DARCY_EQUATION_TYPE, &
    & OC_PROBLEM_STANDARD_DARCY_SUBTYPE],problem,err)
  !Finish the creation of a problem.
  CALL OC_Problem_CreateFinish(problem,err)
  !Start the creation of the problem control loop
  CALL OC_Problem_ControlLoopCreateStart(problem,err)
  !Finish creating the problem control loop
  CALL OC_Problem_ControlLoopCreateFinish(problem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! SOLVER
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of the problem solvers
  CALL OC_Solver_Initialise(linearSolverDarcy,err)
  CALL OC_Problem_SolversCreateStart(problem,err)
  !Get the Darcy solver
  CALL OC_Problem_SolverGet(problem,OC_CONTROL_LOOP_NODE,1,linearSolverDarcy,err)
  !Set the output type
  CALL OC_Solver_OutputTypeSet(linearSolverDarcy,linearSolverDarcyOutputType,err)
  !Set the solver settings
  IF(linearSolverDarcyDirectFlag) THEN
    CALL OC_Solver_LinearTypeSet(linearSolverDarcy,OC_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,err)
    CALL OC_Solver_LibraryTypeSet(linearSolverDarcy,OC_SOLVER_MUMPS_LIBRARY,err)
  ELSE
    CALL OC_Solver_LinearTypeSet(linearSolverDarcy,OC_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,err)
    CALL OC_Solver_LinearIterativeMaximumIterationsSet(linearSolverDarcy,maximumIterations,err)
    CALL OC_Solver_LinearIterativeDivergenceToleranceSet(linearSolverDarcy,divergenceTolerance,err)
    CALL OC_Solver_LinearIterativeRelativeToleranceSet(linearSolverDarcy,relativeTolerance,err)
    CALL OC_Solver_LinearIterativeAbsoluteToleranceSet(linearSolverDarcy,absoluteTolerance,err)
    CALL OC_Solver_LinearIterativeGMRESRestartSet(linearSolverDarcy,restartValue,err)
  ENDIF
  !Finish the creation of the problem solver
  CALL OC_Problem_SolversCreateFinish(problem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! SOLVER EQUATIONS
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of the problem solver equations
  CALL OC_Solver_Initialise(linearSolverDarcy,err)
  CALL OC_SolverEquations_Initialise(solverEquationsDarcy,err)

  CALL OC_Problem_SolverEquationsCreateStart(problem,err)
  !Get the Darcy solver equations
  CALL OC_Problem_SolverGet(problem,OC_CONTROL_LOOP_NODE,1,linearSolverDarcy,err)
  CALL OC_Solver_SolverEquationsGet(linearSolverDarcy,solverEquationsDarcy,err)
  !Set the solver equations sparsity
  CALL OC_SolverEquations_SparsityTypeSet(solverEquationsDarcy,OC_SOLVER_SPARSE_MATRICES,err)
  !Add in the equations set
  CALL OC_SolverEquations_EquationsSetAdd(solverEquationsDarcy,equationsSetDarcy,equationsSetIndex,err)
  !Finish the creation of the problem solver equations
  CALL OC_Problem_SolverEquationsCreateFinish(problem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! BOUNDARY CONDITIONS
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of the equations set boundary conditions for Darcy
  CALL OC_BoundaryConditions_Initialise(boundaryConditionsDarcy,err)
  CALL OC_SolverEquations_BoundaryConditionsCreateStart(solverEquationsDarcy,boundaryConditionsDarcy,err)
  !Set fixed wall nodes
  IF(fixedWallNodesFlag) THEN
    DO nodeIdx=1,numberOfFixedWallNodes
      nodeNumber=fixedWallNodes(nodeIdx)
      condition=OC_BOUNDARY_CONDITION_FIXED
      CALL OC_Decomposition_NodeDomainGet(decomposition,nodeNumber,1,boundaryNodeDomain,err)
      IF(boundaryNodeDomain==computationalNodeNumber) THEN
        DO componentIdx=1,numberOfDimensions
          boundaryConditionValue=0.0_OC_RP
          CALL OC_BoundaryConditions_SetNode(boundaryConditionsDarcy,dependentFieldDarcy,OC_FIELD_U_VARIABLE_TYPE,1, &
            & OC_NO_GLOBAL_DERIV,nodeNumber,componentIdx,condition,boundaryConditionValue,err)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
  !Set velocity boundary conditions
  IF(inletWallNodesFlag) THEN
    DO nodeIdx=1,numberOfInletWallNodes
      nodeNumber=inletWallNodes(nodeIdx)
      condition=OC_BOUNDARY_CONDITION_FIXED
      CALL OC_Decomposition_NodeDomainGet(decomposition,nodeNumber,1,boundaryNodeDomain,err)
      IF(boundaryNodeDomain==computationalNodeNumber) THEN
        DO componentIdx=1,numberOfDimensions
          boundaryConditionValue=boundaryConditionValues(componentIdx)
          CALL OC_BoundaryConditions_SetNode(boundaryConditionsDarcy,dependentFieldDarcy,OC_FIELD_U_VARIABLE_TYPE,1, &
            & OC_NO_GLOBAL_DERIV,nodeNumber,componentIdx,condition,boundaryConditionValue,err)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
  !CALL OC_SolverEquations_BoundaryConditionsAnalytic(solverEquationsDarcy,err)
  !Finish the creation of the equations set boundary conditions
  CALL OC_SolverEquations_BoundaryConditionsCreateFinish(solverEquationsDarcy,err)

  !-----------------------------------------------------------------------------------------------------------
  ! SOLVE
  !-----------------------------------------------------------------------------------------------------------

  INQUIRE(file="./results", exist=directoryExists)
  IF(.NOT.directoryExists) CALL EXECUTE_COMMAND_LINE("mkdir ./output")
  
  !Solve the problem
  WRITE(*,'(A)') "Solving problem..."
  CALL OC_Problem_Solve(problem,err)
  WRITE(*,'(A)') "Problem solved!"

  !-----------------------------------------------------------------------------------------------------------
  ! OUTPUT
  !-----------------------------------------------------------------------------------------------------------

  !Output Analytic Analysis
  !CALL OC_AnalyticAnalysis_Output(dependentFieldDarcy,"DarcyAnalytic",err)

  exportFieldIO=.FALSE.
  IF(exportFieldIO) THEN
     WRITE(*,'(A)') "Exporting fields..."
     CALL OC_Fields_Initialise(fields,err)
     CALL OC_Fields_Create(region,fields,err)
     CALL OC_Fields_NodesExport(fields,"darcy_static","FORTRAN",err)
     CALL OC_Fields_ElementsExport(fields,"darcy_static","FORTRAN",err)
     CALL OC_Fields_Finalise(fields,err)
     WRITE(*,'(A)') "Field exported!"
  ENDIF

  !Destroy the context
  CALL OC_Context_Destroy(context,err)
  !Finialise OpenCMISS
  CALL OC_Finalise(err)
  
  WRITE(*,'(A)') "Program successfully completed."

END PROGRAM DarcyStaticExample
