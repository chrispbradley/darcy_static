import os

# Intialise OpenCMISS-Iron
from opencmiss.iron import iron

#parameters.parse()

#-----------------------------------------------------------------------------------------------------------
# SET PROBLEM PARAMETERS
#-----------------------------------------------------------------------------------------------------------

# Set problem parameters
height = 1.0
width = 1.0
length = 1.0

porosity = 0.3
perm_over_vis = 0.8
initial_conc = 0.0
screen_output_freq = 2 #how many time steps between outputs to screen

(contextUserNumber,
 coordinateSystemUserNumber,
 regionUserNumber,
 basisUserNumber,
 generatedMeshUserNumber,
 meshUserNumber,
 decompositionUserNumber,
 decomposerUserNumber,
 geometricFieldUserNumber,
 equationsSetFieldUserNumber,
 dependentFieldUserNumber,
 materialFieldUserNumber,
 equationsSetUserNumber,
 problemUserNumber) = range(1,15)

numberGlobalXElements = 5
numberGlobalYElements = 5
numberGlobalZElements = 5

context = iron.Context()
context.Create(contextUserNumber)

worldRegion = iron.Region()
context.WorldRegionGet(worldRegion)

#-----------------------------------------------------------------------------------------------------------
# DIAGNOSTICS AND COMPUTATIONAL NODE INFORMATION
#-----------------------------------------------------------------------------------------------------------

computationEnvironment = iron.ComputationEnvironment()
context.ComputationEnvironmentGet(computationEnvironment)

worldWorkGroup = iron.WorkGroup()
computationEnvironment.WorldWorkGroupGet(worldWorkGroup)
numberOfComputationalNodes = worldWorkGroup.NumberOfGroupNodesGet()
computationalNodeNumber = worldWorkGroup.GroupNodeNumberGet()

#-----------------------------------------------------------------------------------------------------------
#COORDINATE SYSTEM
#-----------------------------------------------------------------------------------------------------------

# Create a RC coordinate system
coordinateSystem = iron.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber,context)
coordinateSystem.dimension = 3
coordinateSystem.CreateFinish()

#-----------------------------------------------------------------------------------------------------------
#REGION
#-----------------------------------------------------------------------------------------------------------

# Create a region
region = iron.Region()
region.CreateStart(regionUserNumber,worldRegion)
region.label = "DarcyRegion"
region.coordinateSystem = coordinateSystem
region.CreateFinish()

#-----------------------------------------------------------------------------------------------------------
#BASIS
#-----------------------------------------------------------------------------------------------------------

# Create a tri-linear lagrange basis
basis = iron.Basis()
basis.CreateStart(basisUserNumber,context)
basis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
basis.numberOfXi = 3
basis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*3
basis.quadratureNumberOfGaussXi = [3]*3
basis.quadratureLocalFaceGaussEvaluate = True
basis.CreateFinish()

#-----------------------------------------------------------------------------------------------------------
#MESH
#-----------------------------------------------------------------------------------------------------------

# Create a generated mesh
generatedMesh = iron.GeneratedMesh()
generatedMesh.CreateStart(generatedMeshUserNumber,region)
generatedMesh.type = iron.GeneratedMeshTypes.REGULAR
generatedMesh.basis = [basis]
generatedMesh.extent = [width,height,length]
generatedMesh.numberOfElements = [numberGlobalXElements,numberGlobalYElements,numberGlobalZElements]

mesh = iron.Mesh()
generatedMesh.CreateFinish(meshUserNumber,mesh)

numberOfElements = mesh.numberOfElements
print("number of elements: " + str(numberOfElements))

#-----------------------------------------------------------------------------------------------------------
#MESH DECOMPOSITION
#-----------------------------------------------------------------------------------------------------------

# Create a decomposition for the mesh
decomposition = iron.Decomposition()
decomposition.CreateStart(decompositionUserNumber,mesh)
decomposition.CreateFinish()

#-----------------------------------------------------------------------------------------------------------
#DECOMPOSER
#-----------------------------------------------------------------------------------------------------------

decomposer = iron.Decomposer()
decomposer.CreateStart(decomposerUserNumber,worldRegion,worldWorkGroup)
decompositionIndex = decomposer.DecompositionAdd(decomposition)
decomposer.CreateFinish()

#-----------------------------------------------------------------------------------------------------------
#GEOMETRIC FIELD
#-----------------------------------------------------------------------------------------------------------

# Create a field for the geometry
geometricField = iron.Field()
geometricField.CreateStart(geometricFieldUserNumber,region)
geometricField.decomposition = decomposition
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,1)
geometricField.CreateFinish()

# Set geometry from the generated mesh
generatedMesh.GeometricParametersCalculate(geometricField)

#-----------------------------------------------------------------------------------------------------------
#EQUATION SETS
#-----------------------------------------------------------------------------------------------------------

# Create standard Diffusion equations set
equationsSetField = iron.Field()
equationsSet = iron.EquationsSet()
equationsSetSpecification = [iron.EquationsSetClasses.FLUID_MECHANICS,
        iron.EquationsSetTypes.DARCY_EQUATION,
        iron.EquationsSetSubtypes.STANDARD_DARCY]
equationsSet.CreateStart(equationsSetUserNumber,region,geometricField,
        equationsSetSpecification,equationsSetFieldUserNumber,equationsSetField)
equationsSet.CreateFinish()

#-----------------------------------------------------------------------------------------------------------
#DEPENDENT FIELD
#-----------------------------------------------------------------------------------------------------------

# Create dependent field
dependentField = iron.Field()
equationsSet.DependentCreateStart(dependentFieldUserNumber,dependentField)
dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,1)
dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,4,1)
equationsSet.DependentCreateFinish()

#-----------------------------------------------------------------------------------------------------------
#MATERIAL FIELD
#-----------------------------------------------------------------------------------------------------------

# Create material field
materialField = iron.Field()
equationsSet.MaterialsCreateStart(materialFieldUserNumber,materialField)
materialField.VariableLabelSet(iron.FieldVariableTypes.U, "Material")
equationsSet.MaterialsCreateFinish()

## I believe this will change the diffusion coeff
# k and mu?
materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,porosity)
materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,perm_over_vis)
#materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3,diff_coeff)

# Initialise dependent field
dependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,initial_conc)

#-----------------------------------------------------------------------------------------------------------
# EQUATIONS
#-----------------------------------------------------------------------------------------------------------

# Create equations
equations = iron.Equations()
equationsSet.EquationsCreateStart(equations)
equations.sparsityType = iron.EquationsSparsityTypes.SPARSE
equations.outputType = iron.EquationsOutputTypes.NONE
equationsSet.EquationsCreateFinish()

#-----------------------------------------------------------------------------------------------------------
#PROBLEM
#-----------------------------------------------------------------------------------------------------------

# Create Diffusion equation problem
problem = iron.Problem()
problemSpecification = [iron.ProblemClasses.FLUID_MECHANICS,
        iron.ProblemTypes.DARCY_EQUATION,
        iron.ProblemSubtypes.STANDARD_DARCY]
problem.CreateStart(problemUserNumber,context,problemSpecification)
problem.CreateFinish()

# Create control loops
problem.ControlLoopCreateStart()
problem.ControlLoopCreateFinish()

#-----------------------------------------------------------------------------------------------------------
#SOLVER
#-----------------------------------------------------------------------------------------------------------

# Create problem solver
solver = iron.Solver()
problem.SolversCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,solver)
solver.outputType = iron.SolverOutputTypes.SOLVER
solver.linearType = iron.LinearSolverTypes.ITERATIVE
solver.linearIterativeAbsoluteTolerance = 1.0E-12
solver.linearIterativeRelativeTolerance = 1.0E-12
problem.SolversCreateFinish()

#-----------------------------------------------------------------------------------------------------------
#SOLVER EQUATIONS
#-----------------------------------------------------------------------------------------------------------

## Create solver equations and add equations set to solver equations
solver = iron.Solver()
solverEquations = iron.SolverEquations()
problem.SolverEquationsCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,solver)
solver.SolverEquationsGet(solverEquations)
solverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()

#-----------------------------------------------------------------------------------------------------------
#BOUNDARY CONDITIONS
#-----------------------------------------------------------------------------------------------------------

## Create boundary conditions and set first and last nodes to 0.0 and 1.0
boundaryConditions = iron.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)
firstNodeNumber=1
nodes = iron.Nodes()
region.NodesGet(nodes)
lastNodeNumber = nodes.numberOfNodes
firstNodeDomain = decomposition.NodeDomainGet(firstNodeNumber,1)
lastNodeDomain = decomposition.NodeDomainGet(lastNodeNumber,1)

for i in range(1,37):
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,i,3,iron.BoundaryConditionsTypes.FIXED,1.0)
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,i,1,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,i,2,iron.BoundaryConditionsTypes.FIXED,0.0)
for i in range(181,217):
#if lastNodeDomain == computationalNodeNumber:
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,i,3,iron.BoundaryConditionsTypes.FIXED,1.0)
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,i,1,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,i,2,iron.BoundaryConditionsTypes.FIXED,0.0)
for i in range(145,151):
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,i,2,iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
for i in range(109,115):
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,i,2,iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
for i in range(73,79):
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,i,2,iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
for i in range(37,43):
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,i,2,iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
for i in range(1,7):
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,i,2,iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
for i in range(67,73):
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,i,2,iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
for i in range(31,37):
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,i,2,iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
for i in range(103,109):
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,i,2,iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
for i in range(139,145):
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,i,2,iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
for i in range(175,181):
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,i,2,iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
for i in range(211,217):
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,i,2,iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
mycounter=1
mycounter2=6
for i in range(0,36):
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,mycounter,1,iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,mycounter2,1,iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
    mycounter=mycounter+6
    mycounter2=mycounter2+6
#	1 7 13 19 25 31 37 43 49 55 61 67

solverEquations.BoundaryConditionsCreateFinish()

#-----------------------------------------------------------------------------------------------------------
#SOLVE
#-----------------------------------------------------------------------------------------------------------

problem.Solve()

#-----------------------------------------------------------------------------------------------------------
#OUTPUT
#-----------------------------------------------------------------------------------------------------------

# Export results
#baseName = "Diffusion"
#dataFormat = "PLAIN_TEXT"
#fml = iron.FieldMLIO()
#fml.OutputCreate(mesh, "", baseName, dataFormat)
#fml.OutputAddFieldNoType(baseName+".geometric", dataFormat, geometricField,
#    iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
#fml.OutputAddFieldNoType(baseName+".phi", dataFormat, dependentField,
#    iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
#fml.OutputWrite("DiffusionExample.xml")
#fml.Finalise()

# Ensure output directories exist
if not os.path.exists('./output'):
    os.makedirs('./output')

# Export results
fields = iron.Fields()
fields.CreateRegion(region)
fields.NodesExport("output/StaticDarcy","FORTRAN")
fields.ElementsExport("output/StaticDarcy","FORTRAN")
fields.Finalise()

context.Destroy(context)

Finalise()
