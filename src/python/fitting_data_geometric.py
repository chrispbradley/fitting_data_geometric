#!/usr/bin/env python

#> \file
#> \author Chris Bradley
#> \brief This is an example to use linear fitting to fit a cube mesh surface to a sphere.
#>

import sys, os
import exfile
import numpy
import math
import random

# Intialise OpenCMISS
from opencmiss.opencmiss import OpenCMISS_Python as oc

# defining the output file to be written in the ExDataFile
def writeExdataFile(filename,dataPointLocations,dataErrorVector,dataErrorDistance,offset):
    "Writes data points to an exdata file"

    numberOfDimensions = dataPointLocations[1].shape[0]
    try:
        f = open(filename,"w")
        if numberOfDimensions == 1:
            header = '''Group name: DataPoints
 #Fields=3
 1) data_coordinates, coordinate, rectangular cartesian, #Components='''+str(numberOfDimensions)+'''
  x.  Value index=1, #Derivatives=0, #Versions=1
 2) data_error, field, rectangular cartesian, #Components='''+str(numberOfDimensions)+'''
  x.  Value index=2, #Derivatives=0, #Versions=1
 3) data_distance, field, real, #Components=1
  1.  Value index=3, #Derivatives=0, #Versions=1
'''
        elif numberOfDimensions == 2:
            header = '''Group name: DataPoints
 #Fields=3
 1) data_coordinates, coordinate, rectangular cartesian, #Components='''+str(numberOfDimensions)+'''
  x.  Value index=1, #Derivatives=0, #Versions=1
  y.  Value index=2, #Derivatives=0, #Versions=1
 2) data_error, field, rectangular cartesian, #Components='''+str(numberOfDimensions)+'''
  x.  Value index=3, #Derivatives=0, #Versions=1
  y.  Value index=4, #Derivatives=0, #Versions=1
 3) data_distance, field, real, #Components=1
  1.  Value index=5, #Derivatives=0, #Versions=1
'''
        elif numberOfDimensions == 3:
             header = '''Group name: DataPoints
 #Fields=3
 1) data_coordinates, coordinate, rectangular cartesian, #Components='''+str(numberOfDimensions)+'''
  x.  Value index=1, #Derivatives=0, #Versions=1
  y.  Value index=2, #Derivatives=0, #Versions=1
  x.  Value index=3, #Derivatives=0, #Versions=1
 2) data_error, field, rectangular cartesian, #Components='''+str(numberOfDimensions)+'''
  x.  Value index=4, #Derivatives=0, #Versions=1
  y.  Value index=5, #Derivatives=0, #Versions=1
  z.  Value index=6, #Derivatives=0, #Versions=1
 3) data_distance, field, real, #Components=1
  1.  Value index=7, #Derivatives=0, #Versions=1
'''
        f.write(header)

        numberOfDataPoints = len(dataPointLocations)
        for i in range(numberOfDataPoints):
            line = " Node: " + str(offset+i+1) + '\n'
            f.write(line)
            for j in range (numberOfDimensions):
                line = ' ' + str(dataPointLocations[i,j]) + '\t'
                f.write(line)
            line = '\n'
            f.write(line)
            for j in range (numberOfDimensions):
                line = ' ' + str(dataErrorVector[i,j]) + '\t'
                f.write(line)
            line = '\n'
            f.write(line)
            line = ' ' + str(dataErrorDistance[i])
            f.write(line)
            line = '\n'
            f.write(line)
        f.close()
            
    except IOError:
        print ('Could not open file: ' + filename)

#=================================================================
# Control Panel
#=================================================================

# Set data point resolution (will be randomly placed on surface of a sphere)
numberOfDataPoints = 30
# Set Sobolev smoothing parameters
tau = 5.0
kappa = 1.0
# iteratively fit the cube to sphere- default 1 for automated testing
numberOfIterations = 3

# Override with command line arguments if need be
if len(sys.argv) > 1:
    if len(sys.argv) > 5:
        sys.exit('Error: too many arguments- currently only accepting 4 options: numberOfDataPoints tau kappa numberOfIterations')
    numberOfDataPoints = int(sys.argv[1])
    if len(sys.argv) > 2:
        tau = float(sys.argv[2])
    if len(sys.argv) > 3:
        kappa = float(sys.argv[3])
    if len(sys.argv) > 4:
        numberOfIterations = int(sys.argv[4])

# Set cube dimensions
CUBE_SIZE = 1.00

NUMBER_OF_GAUSS_XI = 3
ZERO_TOLERANCE = 0.00001
DATA_OFFSET=1000

#=================================================================

(CONTEXT_USER_NUMBER,
 COORDINATE_SYSTEM_USER_NUMBER,
 REGION_USER_NUMBER,
 BASIS_USER_NUMBER,
 MESH_USER_NUMBER,
 DECOMPOSITION_USER_NUMBER,
 DECOMPOSER_USER_NUMBER,
 GEOMETRIC_FIELD_USER_NUMBER,
 EQUATIONS_SET_FIELD_USER_NUMBER,
 DEPENDENT_FIELD_USER_NUMBER,
 INDEPENDENT_FIELD_USER_NUMBER,
 DATA_POINT_USER_NUMBER,
 MATERIALS_FIELD_USER_NUMBER,
 DATA_PROJECTION_USER_NUMBER,
 EQUATIONS_SET_USER_NUMBER,
 PROBLEM_USER_NUMBER) = range(1,17)

context = oc.Context()
context.Create(CONTEXT_USER_NUMBER)

worldRegion = oc.Region()
context.WorldRegionGet(worldRegion)

# Get the computational nodes information
computationEnvironment = oc.ComputationEnvironment()
context.ComputationEnvironmentGet(computationEnvironment)

worldWorkGroup = oc.WorkGroup()
computationEnvironment.WorldWorkGroupGet(worldWorkGroup)
numberOfComputationalNodes = worldWorkGroup.NumberOfGroupNodesGet()
computationalNodeNumber = worldWorkGroup.GroupNodeNumberGet()

# Create a RC coordinate system
coordinateSystem = oc.CoordinateSystem()
coordinateSystem.CreateStart(COORDINATE_SYSTEM_USER_NUMBER,context)
coordinateSystem.DimensionSet(3)
coordinateSystem.CreateFinish()

# Create a region
region = oc.Region()
region.CreateStart(REGION_USER_NUMBER,worldRegion)
region.LabelSet("FittingRegion")
region.CoordinateSystemSet(coordinateSystem)
region.CreateFinish()

#=================================================================
# Mesh
#=================================================================

# Create a tricubic Hermite basis
basis = oc.Basis()
basis.CreateStart(BASIS_USER_NUMBER,context)
basis.TypeSet(oc.BasisTypes.LAGRANGE_HERMITE_TP)
basis.NumberOfXiSet(3)
basis.InterpolationXiSet([oc.BasisInterpolationSpecifications.CUBIC_HERMITE]*3)
basis.QuadratureNumberOfGaussXiSet([NUMBER_OF_GAUSS_XI]*3)
basis.CreateFinish()

# Define nodes for the mesh
nodes = oc.Nodes()
nodes.CreateStart(region,8)
nodes.CreateFinish()

# Create the mesh
mesh = oc.Mesh()
mesh.CreateStart(MESH_USER_NUMBER,region,3)
mesh.NumberOfComponentsSet(1)
mesh.NumberOfElementsSet(1)
elements = oc.MeshElements()
elements.CreateStart(mesh, 1, basis)
elements.NodesSet(1,[1,2,3,4,5,6,7,8])
elements.CreateFinish()
mesh.CreateFinish()

# Create a decomposition for the mesh
decomposition = oc.Decomposition()
decomposition.CreateStart(DECOMPOSITION_USER_NUMBER,mesh)
decomposition.CalculateFacesSet(True)
decomposition.CreateFinish()

# Decompose 
decomposer = oc.Decomposer()
decomposer.CreateStart(DECOMPOSER_USER_NUMBER,worldRegion,worldWorkGroup)
decompositionIndex = decomposer.DecompositionAdd(decomposition)
decomposer.CreateFinish()

#=================================================================
# Geometric Field
#=================================================================

# Create a field for the geometry
geometricField = oc.Field()
geometricField.CreateStart(GEOMETRIC_FIELD_USER_NUMBER,region)
geometricField.DecompositionSet(decomposition)
geometricField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,1,1)
geometricField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,2,1)
geometricField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,3,1)
geometricField.ScalingTypeSet(oc.FieldScalingTypes.ARITHMETIC_MEAN)
geometricField.CreateFinish()

# Set the geometric field dofs
# Node 1
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 1, 1, -CUBE_SIZE)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 1, 2, -CUBE_SIZE)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 1, 3, -CUBE_SIZE)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 1, 1, 1.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 1, 2, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 1, 3, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 1, 1, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 1, 2, 1.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 1, 3, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 1, 1, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 1, 2, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 1, 3, 1.0)
# Node 2
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 2, 1,  CUBE_SIZE)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 2, 2, -CUBE_SIZE)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 2, 3, -CUBE_SIZE)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 2, 1, 1.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 2, 2, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 2, 3, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 2, 1, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 2, 2, 1.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 2, 3, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 2, 1, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 2, 2, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 2, 3, 1.0)
# Node 3
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 3, 1, -CUBE_SIZE)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 3, 2,  CUBE_SIZE)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 3, 3, -CUBE_SIZE)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 3, 1, 1.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 3, 2, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 3, 3, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 3, 1, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 3, 2, 1.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 3, 3, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 3, 1, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 3, 2, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 3, 3, 1.0)
# Node 4
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 4, 1,  CUBE_SIZE)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 4, 2,  CUBE_SIZE)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 4, 3, -CUBE_SIZE)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 4, 1, 1.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 4, 2, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 4, 3, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 4, 1, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 4, 2, 1.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 4, 3, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 4, 1, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 4, 2, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 4, 3, 1.0)
# Node 5
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 5, 1, -CUBE_SIZE)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 5, 2, -CUBE_SIZE)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 5, 3,  CUBE_SIZE)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 5, 1, 1.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 5, 2, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 5, 3, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 5, 1, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 5, 2, 1.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 5, 3, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 5, 1, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 5, 2, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 5, 3, 1.0)
# Node 6
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 6, 1,  CUBE_SIZE)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 6, 2, -CUBE_SIZE)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 6, 3,  CUBE_SIZE)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 6, 1, 1.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 6, 2, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 6, 3, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 6, 1, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 6, 2, 1.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 6, 3, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 6, 1, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 6, 2, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 6, 3, 1.0)
# Node 7
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 7, 1, -CUBE_SIZE)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 7, 2,  CUBE_SIZE)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 7, 3,  CUBE_SIZE)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 7, 1, 1.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 7, 2, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 7, 3, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 7, 1, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 7, 2, 1.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 7, 3, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 7, 1, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 7, 2, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 7, 3, 1.0)
# Node 8
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 8, 1,  CUBE_SIZE)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 8, 2,  CUBE_SIZE)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 8, 3,  CUBE_SIZE)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 8, 1, 1.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 8, 2, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 8, 3, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 8, 1, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 8, 2, 1.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 8, 3, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 8, 1, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 8, 2, 0.0)
geometricField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                      1, oc.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 8, 3, 1.0)

geometricField.ParameterSetUpdateStart(oc.FieldVariableTypes.U, oc.FieldParameterSetTypes.VALUES)
geometricField.ParameterSetUpdateFinish(oc.FieldVariableTypes.U, oc.FieldParameterSetTypes.VALUES)

#=================================================================
# Data Points
#=================================================================

# Create the data points
dataPoints = oc.DataPoints()
dataPoints.CreateStart(DATA_POINT_USER_NUMBER,region,numberOfDataPoints)

localNumberOfDataPoints = 0
dataPointLocations = numpy.zeros((numberOfDataPoints,3))
print("Number of data points: " + str(numberOfDataPoints))

# Calculate data point locations of points on a sphere
random.seed(1)
for i in range(numberOfDataPoints):
    x = 2.0*CUBE_SIZE*(random.uniform(0.0,1.0)-0.5)
    y = 2.0*CUBE_SIZE*(random.uniform(0.0,1.0)-0.5)
    r2 = 3.0*CUBE_SIZE*CUBE_SIZE
    if (r2-x*x-y*y)>ZERO_TOLERANCE :
        z = math.sqrt(r2-x*x-y*y)
    else:
        z = CUBE_SIZE
    dataPointLocations[i,:] = [x,y,z]

# Set up CMISS data points with geometric values
for dataPoint in range(numberOfDataPoints):
    dataPointId = dataPoint + 1
    dataList = dataPointLocations[dataPoint,:]
    dataPoints.PositionSet(dataPointId,dataList)
    dataLabel = str(dataPoint + DATA_OFFSET)
    dataPoints.LabelSet(dataPointId,dataLabel)

dataPoints.CreateFinish()

#=================================================================
# Data Projection on Geometric Field
#=================================================================

print("Projecting data points onto geometric field")
# Set up data projection
dataProjection = oc.DataProjection()
dataProjection.CreateStart(DATA_PROJECTION_USER_NUMBER,dataPoints,geometricField,oc.FieldVariableTypes.U)
dataProjection.ProjectionTypeSet(oc.DataProjectionProjectionTypes.BOUNDARY_FACES)
dataProjection.ProjectionCandidateFacesSet([1],[oc.ElementNormalXiDirections.PLUS_XI3])
dataProjection.CreateFinish()

#dataProjection.ResultElementNumberSet(1,1)
#dataProjection.ResultXiSet(1,[0.1,0.1])

# Evaluate data projection based on geometric field
dataProjection.DataPointsProjectionEvaluate(oc.FieldParameterSetTypes.VALUES)
# Create mesh topology for data projection
mesh.TopologyDataPointsCalculateProjection(dataProjection)
# Create decomposition data projection
decomposition.DataProjectionCalculate()

dataProjection.ResultAnalysisOutput("")

dataErrorVector = numpy.zeros((numberOfDataPoints,3))
dataErrorDistance = numpy.zeros(numberOfDataPoints)
elementIdx=1
numberOfProjectedDataPoints = decomposition.NumberOfElementDataPointsGet(elementIdx)
for dataPointIdx in range(1,numberOfProjectedDataPoints+1):
    dataPointNumber = decomposition.ElementDataPointUserNumberGet(elementIdx,dataPointIdx)
    errorVector = dataProjection.ResultProjectionVectorGet(dataPointNumber,3)
    dataErrorVector[dataPointNumber-1,0]=errorVector[0]
    dataErrorVector[dataPointNumber-1,1]=errorVector[1]
    dataErrorVector[dataPointNumber-1,2]=errorVector[2]
    errorDistance = dataProjection.ResultDistanceGet(dataPointNumber)
    dataErrorDistance[dataPointNumber-1]=errorDistance
 
# write data points to exdata file for CMGUI
writeExdataFile("DataPoints.part"+str(computationalNodeNumber)+".exdata",dataPointLocations,dataErrorVector,
                dataErrorDistance,DATA_OFFSET)
print("Projection complete")

#=================================================================
# Equations Set
#=================================================================

# Create vector fitting equations set
equationsSetField = oc.Field()
equationsSet = oc.EquationsSet()
equationsSetSpecification = [oc.EquationsSetClasses.FITTING,
                             oc.EquationsSetTypes.DATA_FITTING_EQUATION,
                             oc.EquationsSetSubtypes.GENERALISED_DATA_FITTING,
                             oc.EquationsSetFittingSmoothingTypes.SOBOLEV_VALUE]
equationsSet.CreateStart(EQUATIONS_SET_USER_NUMBER,region,geometricField,
        equationsSetSpecification,EQUATIONS_SET_FIELD_USER_NUMBER,equationsSetField)
equationsSet.CreateFinish()

#=================================================================
# Dependent Field
#=================================================================

# Create dependent field (will be deformed fitted values based on data point locations)
dependentField = oc.Field()
equationsSet.DependentCreateStart(DEPENDENT_FIELD_USER_NUMBER,dependentField)
dependentField.VariableLabelSet(oc.FieldVariableTypes.U,"Dependent")
dependentField.NumberOfComponentsSet(oc.FieldVariableTypes.U,3)
dependentField.ScalingTypeSet(oc.FieldScalingTypes.ARITHMETIC_MEAN)
equationsSet.DependentCreateFinish()

# Initialise dependent field to undeformed geometric field
for component in range (1,4):
    geometricField.ParametersToFieldParametersComponentCopy(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                                            component, dependentField, oc.FieldVariableTypes.U,
                                                            oc.FieldParameterSetTypes.VALUES, component)

#=================================================================
# Independent Field
#=================================================================

# Create data point field (independent field, with vector values stored at the data points)
independentField = oc.Field()
equationsSet.IndependentCreateStart(INDEPENDENT_FIELD_USER_NUMBER,independentField)
independentField.VariableLabelSet(oc.FieldVariableTypes.U,"DataPointVector")
independentField.VariableLabelSet(oc.FieldVariableTypes.V,"DataPointWeight")
independentField.NumberOfComponentsSet(oc.FieldVariableTypes.U,3)
independentField.NumberOfComponentsSet(oc.FieldVariableTypes.V,3)
independentField.DataProjectionSet(dataProjection)
equationsSet.IndependentCreateFinish()

# loop over each element's data points and set independent field values to data point locations on surface of the sphere
elementDomain = decomposition.ElementDomainGet(1)
if (elementDomain == computationalNodeNumber):
    numberOfProjectedDataPoints = decomposition.NumberOfElementDataPointsGet(1)
    for dataPoint in range(numberOfProjectedDataPoints):
        dataPointId = dataPoint + 1
        dataPointNumber = decomposition.ElementDataPointUserNumberGet(1,dataPointId)
        dataList = dataPoints.PositionGet(dataPointNumber,3)        
        x = dataList[0]
        y = dataList[1]
        z = dataList[2]
        # set data point field values
        independentField.ParameterSetUpdateElementDataPointDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES, 
                                                              1,dataPointId,1,x)
        independentField.ParameterSetUpdateElementDataPointDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES, 
                                                              1,dataPointId,2,y)
        independentField.ParameterSetUpdateElementDataPointDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                                              1,dataPointId,3,z)

#=================================================================
# Material Field
#=================================================================

# Create material field (Sobolev parameters)
materialField = oc.Field()
equationsSet.MaterialsCreateStart(MATERIALS_FIELD_USER_NUMBER,materialField)
materialField.VariableLabelSet(oc.FieldVariableTypes.U,"SmoothingParameters")
equationsSet.MaterialsCreateFinish()

# Set kappa and tau - Sobolev smoothing parameters
materialField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,1,tau)
materialField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,2,kappa)

#=================================================================
# Equations
#=================================================================

# Create equations
equations = oc.Equations()
equationsSet.EquationsCreateStart(equations)
equations.sparsityType = oc.EquationsSparsityTypes.SPARSE
equations.outputType = oc.EquationsOutputTypes.NONE
equations.outputType = oc.EquationsOutputTypes.MATRIX
equationsSet.EquationsCreateFinish()

#=================================================================
# Problem setup
#=================================================================

# Create fitting problem
problem = oc.Problem()
problemSpecification = [oc.ProblemClasses.FITTING,
                        oc.ProblemTypes.FITTING,
                        oc.ProblemSubtypes.STATIC_LINEAR_FITTING]
problem.CreateStart(PROBLEM_USER_NUMBER,context,problemSpecification)
problem.CreateFinish()

# Create control loops
problem.ControlLoopCreateStart()
problem.ControlLoopCreateFinish()

# Create problem solver
solver = oc.Solver()
problem.SolversCreateStart()
problem.SolverGet([oc.ControlLoopIdentifiers.NODE],1,solver)
#solver.OutputTypeSet(oc.SolverOutputTypes.NONE)
solver.OutputTypeSet(oc.SolverOutputTypes.MATRIX)
solver.LinearTypeSet(oc.LinearSolverTypes.ITERATIVE)
#solver.LibraryTypeSet(oc.SolverLibraries.UMFPACK) # UMFPACK/SUPERLU
solver.LinearIterativeAbsoluteToleranceSet(1.0E-10)
solver.LinearIterativeRelativeToleranceSet(1.0E-05)
problem.SolversCreateFinish()

# Create solver equations and add equations set to solver equations
solver = oc.Solver()
solverEquations = oc.SolverEquations()
problem.SolverEquationsCreateStart()
problem.SolverGet([oc.ControlLoopIdentifiers.NODE],1,solver)
solver.SolverEquationsGet(solverEquations)
solverEquations.SparsityTypeSet(oc.SolverEquationsSparsityTypes.SPARSE)
equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()

#=================================================================
# Boundary Conditions
#=================================================================

# Create boundary conditions and set first and last nodes to 0.0 and 1.0
boundaryConditions = oc.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)

for nodeIdx in range(1,5):
    for componentIdx in range(1,4):
        for derivativeIdx in range(1,9):
            boundaryConditions.AddNode(dependentField,oc.FieldVariableTypes.U,
                                       1,derivativeIdx,nodeIdx,componentIdx,
                                       oc.BoundaryConditionsTypes.FIXED,0.0)
       
solverEquations.BoundaryConditionsCreateFinish()

# Export undeformed mesh geometry
print("Writing undeformed geometry")
fields = oc.Fields()
fields.CreateRegion(region)
fields.NodesExport("UndeformedGeometry","FORTRAN")
fields.ElementsExport("UndeformedGeometry","FORTRAN")
fields.Finalise()

#=================================================================
# S o l v e    a n d    E x p o r t    D a t a
#=================================================================
derivativeVector=[0.0,0.0,0.0,0.0]
for iteration in range (1,numberOfIterations+1):

    # Solve the problem
    print("Solving fitting problem, iteration: " + str(iteration))
    problem.Solve()
    # Normalise derivatives
    for nodeIdx in range(1,9):
      for derivativeIdx in [oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1,
                            oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2,
                            oc.GlobalDerivativeConstants.GLOBAL_DERIV_S3]:
          length=0.0
          for componentIdx in range(1,4):
              derivativeVector[componentIdx]=dependentField.ParameterSetGetNode(oc.FieldVariableTypes.U,
                                                                                oc.FieldParameterSetTypes.VALUES,
                                                                                1,derivativeIdx,nodeIdx,componentIdx)
              length=length + derivativeVector[componentIdx]*derivativeVector[componentIdx]
          length=math.sqrt(length)
          for componentIdx in range(1,4):
              value=derivativeVector[componentIdx]/length
              dependentField.ParameterSetUpdateNode(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                                    1,derivativeIdx,nodeIdx,componentIdx,value)

    # Copy dependent field to geometric 
    for componentIdx in range(1,4):
        dependentField.ParametersToFieldParametersComponentCopy(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                                                componentIdx,geometricField,oc.FieldVariableTypes.U,
                                                                oc.FieldParameterSetTypes.VALUES,componentIdx)
    # Reproject
    dataProjection.DataPointsProjectionEvaluate(oc.FieldParameterSetTypes.VALUES)
    #dataProjection.ResultAnalysisOutput("")
    rmsError=dataProjection.ResultRMSErrorGet()
    print("RMS error = "+ str(rmsError))
    # Export fields
    print("Writing deformed geometry")
    fields = oc.Fields()
    fields.CreateRegion(region)
    fields.NodesExport("DeformedGeometry" + str(iteration),"FORTRAN")
    fields.ElementsExport("DeformedGeometry" + str(iteration),"FORTRAN")
    fields.Finalise()

# Destroy the context
context.Destroy()
# Finalise OpenCMISS
oc.Finalise()
