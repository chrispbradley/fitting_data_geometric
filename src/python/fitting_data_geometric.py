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
from opencmiss.iron import iron


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
numberOfDataPoints = 50
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
cubeSize = 1.00
sphereRadius = 1.00

# Set origin
origin = [0.0,0.0,0.0]

numberOfGaussXi = 3
zeroTolerance = 0.00001

#=================================================================

(coordinateSystemUserNumber,
    regionUserNumber,
    basisUserNumber,
    generatedMeshUserNumber,
    meshUserNumber,
    decompositionUserNumber,
    geometricFieldUserNumber,
    equationsSetFieldUserNumber,
    dependentFieldUserNumber,
    independentFieldUserNumber,
    dataPointUserNumber,
    dataPointFieldUserNumber,
    materialFieldUserNumber,
    analyticFieldUserNumber,
    dependentDataFieldUserNumber,
    dataProjectionUserNumber,
    equationsSetUserNumber,
    problemUserNumber) = range(1,19)

# Get the computational nodes information
numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
computationalNodeNumber = iron.ComputationalNodeNumberGet()

# Create a RC coordinate system
coordinateSystem = iron.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber)
coordinateSystem.dimension = 3
coordinateSystem.CreateFinish()

# Create a region
region = iron.Region()
region.CreateStart(regionUserNumber,iron.WorldRegion)
region.label = "FittingRegion"
region.coordinateSystem = coordinateSystem
region.CreateFinish()

#=================================================================
# Mesh
#=================================================================

# Create a tricubic Hermite basis
basis = iron.Basis()
basis.CreateStart(basisUserNumber)
basis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
basis.numberOfXi = 3
basis.interpolationXi = [iron.BasisInterpolationSpecifications.CUBIC_HERMITE]*3
basis.quadratureNumberOfGaussXi = [numberOfGaussXi]*3
basis.CreateFinish()

# Define nodes for the mesh
nodes = iron.Nodes()
nodes.CreateStart(region,8)
nodes.CreateFinish()

# Create the mesh
mesh = iron.Mesh()
mesh.CreateStart(meshUserNumber, region, 3)
mesh.NumberOfComponentsSet(1)
mesh.NumberOfElementsSet(1)
elements = iron.MeshElements()
meshComponentNumber = 1
elements.CreateStart(mesh, 1, basis)
elements.NodesSet(1, [1,2,3,4,5,6,7,8])
elements.CreateFinish()
mesh.CreateFinish()

# Create a decomposition for the mesh
decomposition = iron.Decomposition()
decomposition.CreateStart(decompositionUserNumber,mesh)
decomposition.type = iron.DecompositionTypes.CALCULATED
decomposition.numberOfDomains = numberOfComputationalNodes
decomposition.CalculateFacesSet(True)
decomposition.CreateFinish()

#=================================================================
# Geometric Field
#=================================================================

# Create a field for the geometry
geometricField = iron.Field()
geometricField.CreateStart(geometricFieldUserNumber,region)
geometricField.meshDecomposition = decomposition
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,1)
geometricField.ScalingTypeSet(iron.FieldScalingTypes.ARITHMETIC_MEAN)
geometricField.CreateFinish()

# Set the geometric field dofs
# Node 1
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 1, 1, -cubeSize)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 1, 2, -cubeSize)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 1, 3, -cubeSize)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 1, 1, 1.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 1, 2, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 1, 3, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 1, 1, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 1, 2, 1.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 1, 3, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 1, 1, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 1, 2, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 1, 3, 1.0)
# Node 2
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 2, 1,  cubeSize)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 2, 2, -cubeSize)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 2, 3, -cubeSize)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 2, 1, 1.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 2, 2, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 2, 3, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 2, 1, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 2, 2, 1.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 2, 3, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 2, 1, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 2, 2, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 2, 3, 1.0)
# Node 3
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 3, 1, -cubeSize)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 3, 2,  cubeSize)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 3, 3, -cubeSize)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 3, 1, 1.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 3, 2, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 3, 3, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 3, 1, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 3, 2, 1.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 3, 3, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 3, 1, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 3, 2, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 3, 3, 1.0)
# Node 4
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 4, 1,  cubeSize)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 4, 2,  cubeSize)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 4, 3, -cubeSize)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 4, 1, 1.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 4, 2, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 4, 3, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 4, 1, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 4, 2, 1.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 4, 3, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 4, 1, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 4, 2, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 4, 3, 1.0)
# Node 5
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 5, 1, -cubeSize)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 5, 2, -cubeSize)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 5, 3,  cubeSize)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 5, 1, 1.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 5, 2, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 5, 3, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 5, 1, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 5, 2, 1.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 5, 3, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 5, 1, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 5, 2, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 5, 3, 1.0)
# Node 6
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 6, 1,  cubeSize)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 6, 2, -cubeSize)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 6, 3,  cubeSize)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 6, 1, 1.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 6, 2, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 6, 3, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 6, 1, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 6, 2, 1.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 6, 3, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 6, 1, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 6, 2, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 6, 3, 1.0)
# Node 7
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 7, 1, -cubeSize)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 7, 2,  cubeSize)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 7, 3,  cubeSize)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 7, 1, 1.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 7, 2, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 7, 3, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 7, 1, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 7, 2, 1.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 7, 3, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 7, 1, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 7, 2, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 7, 3, 1.0)
# Node 8
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 8, 1,  cubeSize)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 8, 2,  cubeSize)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, 8, 3,  cubeSize)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 8, 1, 1.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 8, 2, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, 8, 3, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 8, 1, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 8, 2, 1.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, 8, 3, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 8, 1, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 8, 2, 0.0)
geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                      1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3, 8, 3, 1.0)

geometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
geometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

#=================================================================
# Data Points
#=================================================================

# Create the data points
dataPoints = iron.DataPoints()
dataPoints.CreateStart(dataPointUserNumber,region,numberOfDataPoints)

localNumberOfDataPoints = 0
dataPointLocations = numpy.zeros((numberOfDataPoints,3))
print("Number of data points: " + str(numberOfDataPoints))

# Calculate data point locations of points on a sphere
random.seed(1)
for i in range(numberOfDataPoints):
    x = 2.0*(random.uniform(0.0,1.0)-0.5)
    y = 2.0*(random.uniform(0.0,1.0)-0.5)
    r2 = 3*cubeSize*cubeSize
    if (r2-x*x-y*y)>zeroTolerance :
        z = math.sqrt(r2-x*x-y*y)
    else:
        z = cubeSize
    dataPointLocations[i,:] = [x,y,z]

# Set up CMISS data points with geometric values
for dataPoint in range(numberOfDataPoints):
    dataPointId = dataPoint + 1
    dataList = dataPointLocations[dataPoint,:]
    dataPoints.PositionSet(dataPointId,dataList)

dataPoints.CreateFinish()

#=================================================================
# Data Projection on Geometric Field
#=================================================================

print("Projecting data points onto geometric field")
# Set up data projection
dataProjection = iron.DataProjection()
dataProjection.CreateStart(dataProjectionUserNumber,dataPoints,geometricField,iron.FieldVariableTypes.U)
dataProjection.projectionType = iron.DataProjectionProjectionTypes.BOUNDARY_FACES
dataProjection.ProjectionCandidateFacesSet([1],[iron.ElementNormalXiDirections.PLUS_XI3])
dataProjection.CreateFinish()

# Evaluate data projection based on geometric field
dataProjection.DataPointsProjectionEvaluate(iron.FieldParameterSetTypes.VALUES)
# Create mesh topology for data projection
mesh.TopologyDataPointsCalculateProjection(dataProjection)
# Create decomposition topology for data projection
decomposition.TopologyDataProjectionCalculate()

dataProjection.ResultAnalysisOutput("")

dataErrorVector = numpy.zeros((numberOfDataPoints,3))
dataErrorDistance = numpy.zeros(numberOfDataPoints)
elementIdx=1
numberOfProjectedDataPoints = decomposition.TopologyNumberOfElementDataPointsGet(elementIdx)
for dataPointIdx in range(1,numberOfProjectedDataPoints+1):
    dataPointNumber = decomposition.TopologyElementDataPointUserNumberGet(elementIdx,dataPointIdx)
    errorVector = dataProjection.ResultProjectionVectorGet(dataPointNumber,3)
    dataErrorVector[dataPointNumber-1,0]=errorVector[0]
    dataErrorVector[dataPointNumber-1,1]=errorVector[1]
    dataErrorVector[dataPointNumber-1,2]=errorVector[2]
    errorDistance = dataProjection.ResultDistanceGet(dataPointNumber)
    dataErrorDistance[dataPointNumber-1]=errorDistance
 
# write data points to exdata file for CMGUI
offset = 0
writeExdataFile("DataPoints.part"+str(computationalNodeNumber)+".exdata",dataPointLocations,dataErrorVector,dataErrorDistance,offset)
print("Projection complete")

#=================================================================
# Equations Set
#=================================================================

# Create vector fitting equations set
equationsSetField = iron.Field()
equationsSet = iron.EquationsSet()
equationsSetSpecification = [iron.EquationsSetClasses.FITTING,
                             iron.EquationsSetTypes.DATA_FITTING_EQUATION,
                             iron.EquationsSetSubtypes.DATA_POINT_FITTING,
                             iron.EquationsSetFittingSmoothingTypes.SOBOLEV_VALUE]
equationsSet.CreateStart(equationsSetUserNumber,region,geometricField,
        equationsSetSpecification, equationsSetFieldUserNumber, equationsSetField)
equationsSet.CreateFinish()

#=================================================================
# Dependent Field
#=================================================================

# Create dependent field (will be deformed fitted values based on data point locations)
dependentField = iron.Field()
equationsSet.DependentCreateStart(dependentFieldUserNumber,dependentField)
dependentField.VariableLabelSet(iron.FieldVariableTypes.U,"Dependent")
dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.U,3)
dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.DELUDELN,3)
dependentField.ScalingTypeSet(iron.FieldScalingTypes.ARITHMETIC_MEAN)
equationsSet.DependentCreateFinish()

# Initialise dependent field to undeformed geometric field
for component in range (1,4):
    geometricField.ParametersToFieldParametersComponentCopy(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                            component, dependentField, iron.FieldVariableTypes.U,
                                                            iron.FieldParameterSetTypes.VALUES, component)

#=================================================================
# Independent Field
#=================================================================

# Create data point field (independent field, with vector values stored at the data points)
independentField = iron.Field()
equationsSet.IndependentCreateStart(independentFieldUserNumber,independentField)
independentField.VariableLabelSet(iron.FieldVariableTypes.U,"DataPointVector")
independentField.VariableLabelSet(iron.FieldVariableTypes.V,"DataPointWeight")
independentField.DataProjectionSet(dataProjection)
equationsSet.IndependentCreateFinish()

# loop over each element's data points and set independent field values to data point locations on surface of the sphere
elementDomain = decomposition.ElementDomainGet(1)
if (elementDomain == computationalNodeNumber):
    numberOfProjectedDataPoints = decomposition.TopologyNumberOfElementDataPointsGet(1)
    for dataPoint in range(numberOfProjectedDataPoints):
        dataPointId = dataPoint + 1
        dataPointNumber = decomposition.TopologyElementDataPointUserNumberGet(1,dataPointId)
        dataList = dataPoints.PositionGet(dataPointNumber,3)        
        x = dataList[0]
        y = dataList[1]
        z = dataList[2]
        # set data point field values
        independentField.ParameterSetUpdateElementDataPointDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, 
                                                              1,dataPointId,1,x)
        independentField.ParameterSetUpdateElementDataPointDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, 
                                                              1,dataPointId,2,y)
        independentField.ParameterSetUpdateElementDataPointDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                              1,dataPointId,3,z)

#=================================================================
# Material Field
#=================================================================

# Create material field (Sobolev parameters)
materialField = iron.Field()
equationsSet.MaterialsCreateStart(materialFieldUserNumber,materialField)
materialField.VariableLabelSet(iron.FieldVariableTypes.U,"SmoothingParameters")
equationsSet.MaterialsCreateFinish()

# Set kappa and tau - Sobolev smoothing parameters
materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,tau)
materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,kappa)

#=================================================================
# Equations
#=================================================================

# Create equations
equations = iron.Equations()
equationsSet.EquationsCreateStart(equations)
equations.sparsityType = iron.EquationsSparsityTypes.SPARSE
equations.outputType = iron.EquationsOutputTypes.NONE
#equations.outputType = iron.EquationsOutputTypes.MATRIX
equationsSet.EquationsCreateFinish()

#=================================================================
# Problem setup
#=================================================================

# Create fitting problem
problem = iron.Problem()
problemSpecification = [iron.ProblemClasses.FITTING,
                        iron.ProblemTypes.DATA_FITTING,
                        iron.ProblemSubtypes.STATIC_FITTING]
problem.CreateStart(problemUserNumber, problemSpecification)
problem.CreateFinish()

# Create control loops
problem.ControlLoopCreateStart()
problem.ControlLoopCreateFinish()

# Create problem solver
solver = iron.Solver()
problem.SolversCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,solver)
solver.outputType = iron.SolverOutputTypes.NONE
#solver.outputType = iron.SolverOutputTypes.MATRIX 
solver.linearType = iron.LinearSolverTypes.ITERATIVE
#solver.LibraryTypeSet(iron.SolverLibraries.UMFPACK) # UMFPACK/SUPERLU
solver.linearIterativeAbsoluteTolerance = 1.0E-10
solver.linearIterativeRelativeTolerance = 1.0E-05
problem.SolversCreateFinish()

# Create solver equations and add equations set to solver equations
solver = iron.Solver()
solverEquations = iron.SolverEquations()
problem.SolverEquationsCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,solver)
solver.SolverEquationsGet(solverEquations)
solverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()

#=================================================================
# Boundary Conditions
#=================================================================

# Create boundary conditions and set first and last nodes to 0.0 and 1.0
boundaryConditions = iron.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)

for nodeIdx in range(1,5):
    for componentIdx in range(1,4):
        for derivativeIdx in range(1,9):
            boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
                                       1,derivativeIdx,nodeIdx,componentIdx,
                                       iron.BoundaryConditionsTypes.FIXED,0.0)
       
solverEquations.BoundaryConditionsCreateFinish()

# Export undeformed mesh geometry
print("Writing undeformed geometry")
fields = iron.Fields()
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
      for derivativeIdx in [iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,
                            iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,
                            iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3]:
          length=0.0
          for componentIdx in range(1,4):
              derivativeVector[componentIdx]=dependentField.ParameterSetGetNode(iron.FieldVariableTypes.U,
                                                                                iron.FieldParameterSetTypes.VALUES,
                                                                                1,derivativeIdx,nodeIdx,componentIdx)
              length=length + derivativeVector[componentIdx]*derivativeVector[componentIdx]
          length=math.sqrt(length)
          for componentIdx in range(1,4):
              value=derivativeVector[componentIdx]/length
              dependentField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,derivativeIdx,nodeIdx,componentIdx,value)

    # Copy dependent field to geometric 
    for componentIdx in range(1,4):
        dependentField.ParametersToFieldParametersComponentCopy(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                                componentIdx,geometricField,iron.FieldVariableTypes.U,
                                                                iron.FieldParameterSetTypes.VALUES,componentIdx)
    # Reproject
    dataProjection.DataPointsProjectionEvaluate(iron.FieldParameterSetTypes.VALUES)
    #dataProjection.ResultAnalysisOutput("")
    rmsError=dataProjection.ResultRMSErrorGet()
    print("RMS error = "+ str(rmsError))
    # Export fields
    print("Writing deformed geometry")
    fields = iron.Fields()
    fields.CreateRegion(region)
    fields.NodesExport("DeformedGeometry" + str(iteration),"FORTRAN")
    fields.ElementsExport("DeformedGeometry" + str(iteration),"FORTRAN")
    fields.Finalise()

iron.Finalise()
