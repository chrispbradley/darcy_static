PROGRAM DarcyStaticExample

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
  REAL(CMISSRP), PARAMETER :: HEIGHT=1.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: WIDTH=1.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: LENGTH=1.0_CMISSRP

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
  INTEGER(CMISSIntg), PARAMETER :: ANALYTIC_FIELD_USER_NUMBER=13
  INTEGER(CMISSIntg), PARAMETER :: EQUATIONS_SET_USER_NUMBER=14
  INTEGER(CMISSIntg), PARAMETER :: PROBLEM_USER_NUMBER=15

  !Program variables
  INTEGER(CMISSIntg) :: numberOfGlobalXElements,numberOfGlobalYElements,numberOfGlobalZElements
  INTEGER(CMISSIntg) :: numberOfDimensions
  INTEGER(CMISSIntg) :: basisType
  INTEGER(CMISSIntg) :: numberOfGaussXi
  INTEGER(CMISSIntg) :: basisXiInterpolation
  INTEGER(CMISSIntg) :: maximumIterations
  INTEGER(CMISSIntg) :: restartValue
  INTEGER(CMISSIntg) :: numberOfFixedWallNodes
  INTEGER(CMISSIntg) :: numberOfInletWallNodes
  INTEGER(CMISSIntg) :: equationsDarcyOutput
  INTEGER(CMISSIntg) :: componentIdx
  INTEGER(CMISSIntg) :: nodeNumber
  INTEGER(CMISSIntg) :: nodeIdx
  INTEGER(CMISSIntg) :: condition
  INTEGER(CMISSIntg) :: linearSolverDarcyOutputType
  INTEGER(CMISSIntg) :: analyticalType
  INTEGER, ALLOCATABLE :: fixedWallNodes(:)
  INTEGER, ALLOCATABLE :: inletWallNodes(:)
  REAL(CMISSRP) :: initialFieldValue(3)
  REAL(CMISSRP) :: boundaryConditionValues(3)
  REAL(CMISSRP) :: divergenceTolerance
  REAL(CMISSRP) :: relativeTolerance
  REAL(CMISSRP) :: absoluteTolerance
  REAL(CMISSRP) :: linesearchAlpha
  REAL(CMISSRP) :: boundaryConditionValue
  REAL(CMISSRP) :: porosityParameter,permOverVisParameter
  LOGICAL :: directoryExists
  LOGICAL :: exportFieldIO
  LOGICAL :: linearSolverDarcyDirectFlag
  LOGICAL :: fixedWallNodesFlag
  LOGICAL :: inletWallNodesFlag

  !CMISS variables
  TYPE(cmfe_BasisType) :: basis
  TYPE(cmfe_BoundaryConditionsType) :: boundaryConditionsDarcy
  TYPE(cmfe_ComputationEnvironmentType) :: computationEnvironment
  TYPE(cmfe_ContextType) :: context
  TYPE(cmfe_ControlLoopType) :: controlLoop
  TYPE(cmfe_CoordinateSystemType) :: coordinateSystem
  TYPE(cmfe_DecomposerType) :: decomposer
  TYPE(cmfe_DecompositionType) :: decomposition
  TYPE(cmfe_EquationsType) :: equationsDarcy
  TYPE(cmfe_EquationsSetType) :: equationsSetDarcy
  TYPE(cmfe_FieldType) :: geometricField,analyticField
  TYPE(cmfe_FieldType) :: equationsSetField
  TYPE(cmfe_FieldType) :: dependentFieldDarcy
  TYPE(cmfe_FieldType) :: materialsFieldDarcy
  TYPE(cmfe_FieldsType) :: fields
  TYPE(cmfe_GeneratedMeshType) :: generatedMesh
  TYPE(cmfe_MeshType) :: mesh
  TYPE(cmfe_ProblemType) :: problem
  TYPE(cmfe_SolverType) :: linearSolverDarcy
  TYPE(cmfe_SolverEquationsType) :: solverEquationsDarcy
  TYPE(cmfe_RegionType) :: region
  TYPE(cmfe_RegionType) :: WorldRegion
  TYPE(cmfe_WorkGroupType) :: worldWorkGroup

  !Generic CMISS variables
  INTEGER(CMISSIntg) :: numberOfComputationalNodes,computationalNodeNumber
  INTEGER(CMISSIntg) :: decompositionIndex,equationsSetIndex,boundaryNodeDomain,err

  !Intialise OpenCMISS
  CALL cmfe_Initialise(err)
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,err)
  !Create a context
  CALL cmfe_Context_Initialise(context,err)
  CALL cmfe_Context_Create(CONTEXT_USER_NUMBER,context,err)
  CALL cmfe_Region_Initialise(worldRegion,err)
  CALL cmfe_Context_WorldRegionGet(context,worldRegion,err)
 
  !Get the computational nodes information
  CALL cmfe_ComputationEnvironment_Initialise(computationEnvironment,err)
  CALL cmfe_Context_ComputationEnvironmentGet(context,computationEnvironment,err)
  
  CALL cmfe_WorkGroup_Initialise(worldWorkGroup,err)
  CALL cmfe_ComputationEnvironment_WorldWorkGroupGet(computationEnvironment,worldWorkGroup,err)
  CALL cmfe_WorkGroup_NumberOfGroupNodesGet(worldWorkGroup,numberOfComputationalNodes,err)
  CALL cmfe_WorkGroup_GroupNodeNumberGet(worldWorkGroup,computationalNodeNumber,err)

  !-----------------------------------------------------------------------------------------------------------
  ! PROBLEM CONTROL PANEL
  !-----------------------------------------------------------------------------------------------------------

  numberOfGlobalXElements=3
  numberOfGlobalYElements=3
  numberOfGlobalZElements=3
  numberOfDimensions=3
  basisType=CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  basisXiInterpolation=CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  !Set initial values
  initialFieldValue(1)=0.0_CMISSRP
  initialFieldValue(2)=0.0_CMISSRP
  initialFieldValue(3)=0.0_CMISSRP
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
    boundaryConditionValues(1)=0.0_CMISSRP
    boundaryConditionValues(2)=1.0_CMISSRP
    boundaryConditionValues(3)=0.0_CMISSRP
  ENDIF
  !Set material parameters
  porosityParameter=0.5_CMISSRP
  permOverVisParameter=0.8_CMISSRP
  !Set number of Gauss points
  numberOfGaussXi=3
  !Set output parameter
  !(NoOutput/ProgressOutput/TimingOutput/SolverOutput/SolverMatrixOutput)
  linearSolverDarcyOutputType=CMFE_SOLVER_NO_OUTPUT
  !(NoOutput/TimingOutput/MatrixOutput/ElementOutput)
  equationsDarcyOutput=CMFE_EQUATIONS_NO_OUTPUT
  !Set solver parameters
  linearSolverDarcyDirectFlag=.FALSE.
  relativeTolerance=1.0E-10_CMISSRP !default: 1.0E-05_CMISSRP
  absoluteTolerance=1.0E-10_CMISSRP !default: 1.0E-10_CMISSRP
  divergenceTolerance=1.0E5_CMISSRP !default: 1.0E5
  maximumIterations=10000_CMISSIntg !default: 100000
  restartValue=3000_CMISSIntg !default: 30
  linesearchAlpha=1.0_CMISSRP

  !IF(numberOfDimensions==2) analyticalType=CMFE_EQUATIONS_SET_DARCY_EQUATION_TWO_DIM_2
  IF(numberOfDimensions==2) analyticalType=CMFE_EQUATIONS_SET_DARCY_EQUATION_TWO_DIM_3
  !IF(numberOfDimensions==3) analyticalType=CMFE_EQUATIONS_SET_DARCY_EQUATION_THREE_DIM_2
  IF(numberOfDimensions==3) analyticalType=CMFE_EQUATIONS_SET_DARCY_EQUATION_THREE_DIM_3

  !-----------------------------------------------------------------------------------------------------------
  ! COORDINATE SYSTEM
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of a new RC coordinate system
  CALL cmfe_CoordinateSystem_Initialise(coordinateSystem,err)
  CALL cmfe_CoordinateSystem_CreateStart(COORDINATE_SYSTEM_USER_NUMBER,context,coordinateSystem,err)
  !Set the coordinate system dimension
  CALL cmfe_CoordinateSystem_DimensionSet(coordinateSystem,numberOfDimensions,err)
  !Finish the creation of the coordinate system
  CALL cmfe_CoordinateSystem_CreateFinish(coordinateSystem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! REGION
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of a new region
  CALL cmfe_Region_Initialise(region,err)
  CALL cmfe_Region_CreateStart(REGION_USER_NUMBER,worldRegion,region,err)
  CALL cmfe_Region_LabelSet(region,"DarcyRegion",err)
  !Set the regions coordinate system as defined above
  CALL cmfe_Region_CoordinateSystemSet(region,coordinateSystem,err)
  !Finish the creation of the region
  CALL cmfe_Region_CreateFinish(region,err)

  !-----------------------------------------------------------------------------------------------------------
  ! BASIS
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of new bases
  CALL cmfe_Basis_Initialise(basis,err)
  CALL cmfe_Basis_CreateStart(BASIS_USER_NUMBER,context,basis,err)
  !Set the basis type (Lagrange/Simplex)
  CALL cmfe_Basis_TypeSet(basis,basisType,err)
  !Set the basis xi number
  CALL cmfe_Basis_NumberOfXiSet(basis,numberOfDimensions,err)
  !Set the basis xi interpolation and number of Gauss points
  IF(numberOfDimensions==2) THEN
    CALL cmfe_Basis_InterpolationXiSet(basis,[basisXiInterpolation,basisXiInterpolation],err)
    CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(basis,[numberOfGaussXi,numberOfGaussXi],err)
  ELSE
    CALL cmfe_Basis_InterpolationXiSet(basis,[basisXiInterpolation,basisXiInterpolation,basisXiInterpolation],err)
    CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(basis,[numberOfGaussXi,numberOfGaussXi,numberOfGaussXi],err)
    CALL cmfe_Basis_QuadratureLocalFaceGaussEvaluateSet(basis,.TRUE.,err)
  ENDIF
  !Finish the creation of the basis
  CALL cmfe_Basis_CreateFinish(basis,err)

  !-----------------------------------------------------------------------------------------------------------
  ! MESH
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of a generated mesh in the region
  CALL cmfe_GeneratedMesh_Initialise(generatedMesh,err)
  CALL cmfe_GeneratedMesh_CreateStart(GENERATED_MESH_USER_NUMBER,region,generatedMesh,err)
  !Set up a regular x*y*z mesh
  CALL cmfe_GeneratedMesh_TypeSet(generatedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,err)
  !Set the default basis
  CALL cmfe_GeneratedMesh_BasisSet(generatedMesh,[basis],err)
  !Define the mesh on the region
  IF(numberOfDimensions==2) THEN
     CALL cmfe_GeneratedMesh_ExtentSet(generatedMesh,[WIDTH,HEIGHT],err)
     CALL cmfe_GeneratedMesh_NumberOfElementsSet(generatedMesh,[numberOfGlobalXElements,numberOfGlobalYElements],err)
  ELSE
     CALL cmfe_GeneratedMesh_ExtentSet(generatedMesh,[WIDTH,HEIGHT,LENGTH],err)
     CALL cmfe_GeneratedMesh_NumberOfElementsSet(generatedMesh,[numberOfGlobalXElements,numberOfGlobalYElements, &
       & numberOfGlobalZElements],err)
  ENDIF
  !Finish the creation of a generated mesh in the region
  CALL cmfe_Mesh_Initialise(mesh,err)
  CALL cmfe_GeneratedMesh_CreateFinish(generatedMesh,MESH_USER_NUMBER,mesh,err)

  !-----------------------------------------------------------------------------------------------------------
  ! MESH DECOMPOSITION
  !-----------------------------------------------------------------------------------------------------------

  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(decomposition,err)
  CALL cmfe_Decomposition_CreateStart(DECOMPOSITION_USER_NUMBER,mesh,decomposition,err)
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
  !Set the field type
  CALL cmfe_Field_TypeSet(geometricField,CMFE_FIELD_GEOMETRIC_TYPE,err)
  !Set the decomposition to use
  CALL cmfe_Field_DecompositionSet(geometricField,decomposition,err)
  !Set the scaling to use
  CALL cmfe_Field_ScalingTypeSet(geometricField,CMFE_FIELD_NO_SCALING,err)
  !Set the mesh component to be used by the field components.
  DO componentIdx=1,numberOfDimensions
     CALL cmfe_Field_ComponentMeshComponentSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,componentIdx, &
       & 1,err)
  ENDDO
  !Finish creating the field
  CALL cmfe_Field_CreateFinish(geometricField,err)
  !Update the geometric field parameters
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(generatedMesh,geometricField,err)

  !-----------------------------------------------------------------------------------------------------------
  ! EQUATIONS SETS
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set for Static Darcy
  CALL cmfe_EquationsSet_Initialise(equationsSetDarcy,err)
  CALL cmfe_Field_Initialise(equationsSetField,err)
  CALL cmfe_EquationsSet_CreateStart(EQUATIONS_SET_USER_NUMBER,region,geometricField,[CMFE_EQUATIONS_SET_FLUID_MECHANICS_CLASS, &
    & CMFE_EQUATIONS_SET_DARCY_EQUATION_TYPE,CMFE_EQUATIONS_SET_STANDARD_DARCY_SUBTYPE],EQUATIONS_SET_FIELD_USER_NUMBER, &
    & equationsSetField,equationsSetDarcy,err)
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(equationsSetDarcy,err)

  !-----------------------------------------------------------------------------------------------------------
  ! DEPENDENT FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set dependent field variables for Static Darcy
  CALL cmfe_Field_Initialise(dependentFieldDarcy,err)
  CALL cmfe_EquationsSet_DependentCreateStart(equationsSetDarcy,DEPENDENT_FIELD_USER_NUMBER,dependentFieldDarcy,err)
  !Set the mesh component to be used by the field components.
  DO componentIdx=1,numberOfDimensions
    CALL cmfe_Field_ComponentMeshComponentSet(dependentFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE,componentIdx,1,err)
    CALL cmfe_Field_ComponentMeshComponentSet(dependentFieldDarcy,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,componentIdx,1,err)
  ENDDO
  componentIdx=numberOfDimensions+1
  CALL cmfe_Field_ComponentMeshComponentSet(dependentFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE,componentIdx,1,err)
  CALL cmfe_Field_ComponentMeshComponentSet(dependentFieldDarcy,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,componentIdx,1,err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(equationsSetDarcy,err)

  !Initialise dependent field (velocity components)
  DO componentIdx=1,numberOfDimensions
    CALL cmfe_Field_ComponentValuesInitialise(dependentFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      & componentIdx,initialFieldValue(componentIdx),err)
  ENDDO
  
  !-----------------------------------------------------------------------------------------------------------
  ! MATERIALS FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set materials field variables for Static Darcy
  CALL cmfe_Field_Initialise(materialsFieldDarcy,err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(equationsSetDarcy,MATERIALS_FIELD_USER_NUMBER,materialsFieldDarcy,err)
  !Finish the equations set materials field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(equationsSetDarcy,err)
  CALL cmfe_Field_ComponentValuesInitialise(materialsFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,porosityParameter,err)
  CALL cmfe_Field_ComponentValuesInitialise(materialsFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 2,permOverVisParameter,err)

  !-----------------------------------------------------------------------------------------------------------
  ! ANALYTIC FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set analytic field variables for static Darcy
  !CALL cmfe_Field_Initialise(analyticField,err)
  !CALL cmfe_EquationsSet_AnalyticCreateStart(equationsSetDarcy,analyticalType,ANALYTIC_FIELD_USER_NUMBER,analyticField,err)
  !Finish the equations set analytic field variables
  !CALL cmfe_EquationsSet_AnalyticCreateFinish(equationsSetDarcy,err)

  !-----------------------------------------------------------------------------------------------------------
  ! EQUATIONS
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set equations
  CALL cmfe_Equations_Initialise(equationsDarcy,err)
  CALL cmfe_EquationsSet_EquationsCreateStart(equationsSetDarcy,equationsDarcy,err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(equationsDarcy,CMFE_EQUATIONS_SPARSE_MATRICES,err)
  !Set the equations set output
  CALL cmfe_Equations_OutputTypeSet(equationsDarcy,equationsDarcyOutput,err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(equationsSetDarcy,err)

  !-----------------------------------------------------------------------------------------------------------
  ! PROBLEM
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of a problem.
  CALL cmfe_Problem_Initialise(problem,err)
  CALL cmfe_ControlLoop_Initialise(controlLoop,err)
  CALL cmfe_Problem_CreateStart(PROBLEM_USER_NUMBER,context,[CMFE_PROBLEM_FLUID_MECHANICS_CLASS,CMFE_PROBLEM_DARCY_EQUATION_TYPE, &
    & CMFE_PROBLEM_STANDARD_DARCY_SUBTYPE],problem,err)
  !Finish the creation of a problem.
  CALL cmfe_Problem_CreateFinish(problem,err)
  !Start the creation of the problem control loop
  CALL cmfe_Problem_ControlLoopCreateStart(problem,err)
  !Finish creating the problem control loop
  CALL cmfe_Problem_ControlLoopCreateFinish(problem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! SOLVER
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of the problem solvers
  CALL cmfe_Solver_Initialise(linearSolverDarcy,err)
  CALL cmfe_Problem_SolversCreateStart(problem,err)
  !Get the Darcy solver
  CALL cmfe_Problem_SolverGet(problem,CMFE_CONTROL_LOOP_NODE,1,linearSolverDarcy,err)
  !Set the output type
  CALL cmfe_Solver_OutputTypeSet(linearSolverDarcy,linearSolverDarcyOutputType,err)
  !Set the solver settings
  IF(linearSolverDarcyDirectFlag) THEN
    CALL cmfe_Solver_LinearTypeSet(linearSolverDarcy,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,err)
    CALL cmfe_Solver_LibraryTypeSet(linearSolverDarcy,CMFE_SOLVER_MUMPS_LIBRARY,err)
  ELSE
    CALL cmfe_Solver_LinearTypeSet(linearSolverDarcy,CMFE_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,err)
    CALL cmfe_Solver_LinearIterativeMaximumIterationsSet(linearSolverDarcy,maximumIterations,err)
    CALL cmfe_Solver_LinearIterativeDivergenceToleranceSet(linearSolverDarcy,divergenceTolerance,err)
    CALL cmfe_Solver_LinearIterativeRelativeToleranceSet(linearSolverDarcy,relativeTolerance,err)
    CALL cmfe_Solver_LinearIterativeAbsoluteToleranceSet(linearSolverDarcy,absoluteTolerance,err)
    CALL cmfe_Solver_LinearIterativeGMRESRestartSet(linearSolverDarcy,restartValue,err)
  ENDIF
  !Finish the creation of the problem solver
  CALL cmfe_Problem_SolversCreateFinish(problem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! SOLVER EQUATIONS
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of the problem solver equations
  CALL cmfe_Solver_Initialise(linearSolverDarcy,err)
  CALL cmfe_SolverEquations_Initialise(solverEquationsDarcy,err)

  CALL cmfe_Problem_SolverEquationsCreateStart(problem,err)
  !Get the Darcy solver equations
  CALL cmfe_Problem_SolverGet(problem,CMFE_CONTROL_LOOP_NODE,1,linearSolverDarcy,err)
  CALL cmfe_Solver_SolverEquationsGet(linearSolverDarcy,solverEquationsDarcy,err)
  !Set the solver equations sparsity
  CALL cmfe_SolverEquations_SparsityTypeSet(solverEquationsDarcy,CMFE_SOLVER_SPARSE_MATRICES,err)
  !Add in the equations set
  CALL cmfe_SolverEquations_EquationsSetAdd(solverEquationsDarcy,equationsSetDarcy,equationsSetIndex,err)
  !Finish the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateFinish(problem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! BOUNDARY CONDITIONS
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of the equations set boundary conditions for Darcy
  CALL cmfe_BoundaryConditions_Initialise(boundaryConditionsDarcy,err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(solverEquationsDarcy,boundaryConditionsDarcy,err)
  !Set fixed wall nodes
  IF(fixedWallNodesFlag) THEN
    DO nodeIdx=1,numberOfFixedWallNodes
      nodeNumber=fixedWallNodes(nodeIdx)
      condition=CMFE_BOUNDARY_CONDITION_FIXED
      CALL cmfe_Decomposition_NodeDomainGet(decomposition,nodeNumber,1,boundaryNodeDomain,err)
      IF(boundaryNodeDomain==computationalNodeNumber) THEN
        DO componentIdx=1,numberOfDimensions
          boundaryConditionValue=0.0_CMISSRP
          CALL cmfe_BoundaryConditions_SetNode(boundaryConditionsDarcy,dependentFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE,1, &
            & CMFE_NO_GLOBAL_DERIV,nodeNumber,componentIdx,condition,boundaryConditionValue,err)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
  !Set velocity boundary conditions
  IF(inletWallNodesFlag) THEN
    DO nodeIdx=1,numberOfInletWallNodes
      nodeNumber=inletWallNodes(nodeIdx)
      condition=CMFE_BOUNDARY_CONDITION_FIXED
      CALL cmfe_Decomposition_NodeDomainGet(decomposition,nodeNumber,1,boundaryNodeDomain,err)
      IF(boundaryNodeDomain==computationalNodeNumber) THEN
        DO componentIdx=1,numberOfDimensions
          boundaryConditionValue=boundaryConditionValues(componentIdx)
          CALL cmfe_BoundaryConditions_SetNode(boundaryConditionsDarcy,dependentFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE,1, &
            & CMFE_NO_GLOBAL_DERIV,nodeNumber,componentIdx,condition,boundaryConditionValue,err)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
  !CALL cmfe_SolverEquations_BoundaryConditionsAnalytic(solverEquationsDarcy,err)
  !Finish the creation of the equations set boundary conditions
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(solverEquationsDarcy,err)

  !-----------------------------------------------------------------------------------------------------------
  ! SOLVE
  !-----------------------------------------------------------------------------------------------------------

  INQUIRE(file="./results", exist=directoryExists)
  IF(.NOT.directoryExists) CALL EXECUTE_COMMAND_LINE("mkdir ./output")
  
  !Solve the problem
  WRITE(*,'(A)') "Solving problem..."
  CALL cmfe_Problem_Solve(problem,err)
  WRITE(*,'(A)') "Problem solved!"

  !-----------------------------------------------------------------------------------------------------------
  ! OUTPUT
  !-----------------------------------------------------------------------------------------------------------

  !Output Analytic Analysis
  !CALL cmfe_AnalyticAnalysis_Output(dependentFieldDarcy,"DarcyAnalytic",err)

  exportFieldIO=.FALSE.
  IF(exportFieldIO) THEN
     WRITE(*,'(A)') "Exporting fields..."
     CALL cmfe_Fields_Initialise(fields,err)
     CALL cmfe_Fields_Create(region,fields,err)
     CALL cmfe_Fields_NodesExport(fields,"darcy_static","FORTRAN",err)
     CALL cmfe_Fields_ElementsExport(fields,"darcy_static","FORTRAN",err)
     CALL cmfe_Fields_Finalise(fields,err)
     WRITE(*,'(A)') "Field exported!"
  ENDIF

  !Destroy the context
  CALL cmfe_Context_Destroy(context,err)
  !Finialise OpenCMISS
  CALL cmfe_Finalise(err)
  
  WRITE(*,'(A)') "Program successfully completed."

END PROGRAM DarcyStaticExample
