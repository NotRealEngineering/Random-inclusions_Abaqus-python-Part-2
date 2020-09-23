#Property of Not Real Engineering 
#Author: Shank S. Kulkarni
#Copyright 2020 Not Real Engineering - All Rights Reserved You may not use, 
#           distribute and modify this code without the written permission 
#           from Not Real Engineering.
############################################################################
##             Creating Random Inclusions                                 ##
############################################################################

from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import random
from array import *
import math
import numpy
import os        # Operating system
import shutil    # copying or moving files

dis=numpy.zeros(1000)

rad=1.0  # radius of inclusion
Max_iterations=11    # Set number of iterations

max_incl = 25      # Set number of inclusions required

for q in range (1,Max_iterations):
    # LET'S CREATE MODEL
    mdb.Model(modelType=STANDARD_EXPLICIT, name='Model-%d' %(q))
    
    # LET'S CREATE PART
    mdb.models['Model-%d' %(q)].ConstrainedSketch(name='__profile__', sheetSize=20.0)
    mdb.models['Model-%d' %(q)].sketches['__profile__'].rectangle(point1=(-10.0, 10.0), 
        point2=(10.0, -10.0))
    mdb.models['Model-%d' %(q)].Part(dimensionality=TWO_D_PLANAR, name='Part-1', type=
        DEFORMABLE_BODY)
    mdb.models['Model-%d' %(q)].parts['Part-1'].BaseShell(sketch=
        mdb.models['Model-%d' %(q)].sketches['__profile__'])
    del mdb.models['Model-%d' %(q)].sketches['__profile__']
    mdb.models['Model-%d' %(q)].ConstrainedSketch(gridSpacing=1.8, name='__profile__', 
        sheetSize=20, transform=
        mdb.models['Model-%d' %(q)].parts['Part-1'].MakeSketchTransform(
        sketchPlane=mdb.models['Model-%d' %(q)].parts['Part-1'].faces[0], 
        sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0.0, 0.0, 0.0)))
    mdb.models['Model-%d' %(q)].parts['Part-1'].projectReferencesOntoSketch(filter=
        COPLANAR_EDGES, sketch=mdb.models['Model-%d' %(q)].sketches['__profile__'])

    num_incl = 0
    x_coordinate = []
    y_coordinate = []

    while (num_incl < max_incl):
        random_x=random.uniform(-8.7, 8.7)  #generate random x_coordinate within RVE
        random_y=random.uniform(-8.7, 8.7)  #generate random y_coordinate within RVE

        isPointIntersecting = False
        # To check if new inclusion intersects with any existing inclusions
        for j in range (0,len(x_coordinate)):
    
    
            dis[j]=sqrt((random_x-x_coordinate[j])**2+(random_y-y_coordinate[j])**2)

                
            if dis[j] < (2.2*rad):

                isPointIntersecting = True
                break

        if (isPointIntersecting == False):
            x_coordinate.append(random_x)
            y_coordinate.append(random_y)
            num_incl = num_incl + 1  # count no of inclusions       


    for i in range(max_incl):    

        mdb.models['Model-%d' %(q)].sketches['__profile__'].CircleByCenterPerimeter(center=(
            x_coordinate[i], y_coordinate[i]), point1=((x_coordinate[i]-rad), y_coordinate[i]))

        mdb.models['Model-%d' %(q)].parts['Part-1'].PartitionFaceBySketch(faces=
            mdb.models['Model-%d' %(q)].parts['Part-1'].faces.findAt(((9.9, 
            9.9, 0.0), (0.0, 0.0, 1.0)), ), sketch=mdb.models['Model-%d' %(q)].sketches['__profile__'])
                    

    # LET'S CREATE MATERIAL-1 (MATRIX POLYMER)
    mdb.models['Model-%d' %(q)].Material(name='Matrix')
    mdb.models['Model-%d' %(q)].materials['Matrix'].Elastic(table=
        ((1e2, 0.47), ))
    
    # LET'S CREATE MATERIAL-2 (ELASTIC INCLUSION)
    mdb.models['Model-%d' %(q)].Material(name='Elastic')
    mdb.models['Model-%d' %(q)].materials['Elastic'].Elastic(table=
        ((1e3, 0.35), ))
        
    # LET'S CREATE SECTIONS    
    mdb.models['Model-%d' %(q)].HomogeneousSolidSection(material='Matrix', name='Matrix', 
        thickness=None)
    mdb.models['Model-%d' %(q)].HomogeneousSolidSection(material='Elastic', name='Inclusion', 
        thickness=None)

    # LET'S ASSIGN SECTIONS
    mdb.models['Model-%d' %(q)].parts['Part-1'].SectionAssignment(offset=0.0, 
        offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
        faces=mdb.models['Model-%d' %(q)].parts['Part-1'].faces.findAt(((9.9, 
        9.9, 0.0), (0.0, 0.0, 1.0)), )), sectionName='Matrix', 
        thicknessAssignment=FROM_SECTION)

    for i in range (num_incl):    
        mdb.models['Model-%d' %(q)].parts['Part-1'].SectionAssignment(offset=0.2, 
            offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
            faces=mdb.models['Model-%d' %(q)].parts['Part-1'].faces.findAt((((x_coordinate[i]), 
            (y_coordinate[i]), 0.0), (0.0, 0.0, 1.0)), )), sectionName='Inclusion', 
            thicknessAssignment=FROM_SECTION)    
  
    # LET'S CREATE INSTANCE
    mdb.models['Model-%d' %(q)].rootAssembly.DatumCsysByDefault(CARTESIAN)
    mdb.models['Model-%d' %(q)].rootAssembly.Instance(dependent=OFF, name='Part-1-1', 
        part=mdb.models['Model-%d' %(q)].parts['Part-1'])
 
    # LET'S CREATE STEP
    mdb.models['Model-%d' %(q)].StaticStep(initialInc=0.01, maxInc=0.1, maxNumInc=10000, 
        minInc=1e-12, name='Step-1', previous='Initial')
        
    # LET'S CREATE BOUNDARY CONDITIONS
    mdb.models['Model-%d' %(q)].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
        distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
        'BC-1', region=Region(
        edges=mdb.models['Model-%d' %(q)].rootAssembly.instances['Part-1-1'].edges.findAt(
        ((-5.0, -10.0, 0.0), ), )), u1=UNSET, u2=0.0, ur3=UNSET)
    mdb.models['Model-%d' %(q)].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
        distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
        'BC-2', region=Region(
        edges=mdb.models['Model-%d' %(q)].rootAssembly.instances['Part-1-1'].edges.findAt(
        ((-10.0, 5.0, 0.0), ), )), u1=0.0, u2=UNSET, ur3=UNSET)
    mdb.models['Model-%d' %(q)].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
        distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
        'BC-3', region=Region(
        edges=mdb.models['Model-%d' %(q)].rootAssembly.instances['Part-1-1'].edges.findAt(
        ((-5.0, 10.0, 0.0), ), )), u1=UNSET, u2=1.0, ur3=UNSET)
    
    # LET'S SEED THE INSTANCE    
    mdb.models['Model-%d' %(q)].rootAssembly.seedPartInstance(deviationFactor=0.1, 
        minSizeFactor=0.1, regions=(
        mdb.models['Model-%d' %(q)].rootAssembly.instances['Part-1-1'], ), size=0.2)
    
    # LET'S SET ELEMENT TYPE
    for i in range (num_incl):
        mdb.models['Model-%d' %(q)].rootAssembly.setElementType(elemTypes=(ElemType(
            elemCode=CPE4R, elemLibrary=STANDARD), ElemType(elemCode=CPE3, 
            elemLibrary=STANDARD)),regions=(
            mdb.models['Model-%d' %(q)].rootAssembly.instances['Part-1-1'].faces.findAt(((
            (x_coordinate[i]), (y_coordinate[i]), 0.0), )), ))        
    
    mdb.models['Model-%d' %(q)].rootAssembly.setElementType(elemTypes=(ElemType(
        elemCode=CPE4R, elemLibrary=STANDARD), ElemType(elemCode=CPE3, 
        elemLibrary=STANDARD)),  regions=(
        mdb.models['Model-%d' %(q)].rootAssembly.instances['Part-1-1'].faces.findAt(((
        9.9, 9.9, 0.0), )), ))
    
    # LET'S GENERATE MESH
    mdb.models['Model-%d' %(q)].rootAssembly.generateMesh(regions=(
        mdb.models['Model-%d' %(q)].rootAssembly.instances['Part-1-1'], ))    

    # LET'S CREATE SETS FOR EDGES
    mdb.models['Model-%d' %(q)].rootAssembly.Set(edges=
        mdb.models['Model-%d' %(q)].rootAssembly.instances['Part-1-1'].edges.findAt(((
        -5.0, -10.0, 0.0), )), name='BottomEdge')
    mdb.models['Model-%d' %(q)].rootAssembly.Set(edges=
        mdb.models['Model-%d' %(q)].rootAssembly.instances['Part-1-1'].edges.findAt(((
        -5.0, 10.0, 0.0), )), name='TopEdge')

    # LET'S CREATE HISTORY OUTPUT REQUESTS
    mdb.models['Model-%d' %(q)].HistoryOutputRequest(createStepName='Step-1', frequency=1
        , name='H-Output-2', rebar=EXCLUDE, region=
        mdb.models['Model-%d' %(q)].rootAssembly.sets['BottomEdge'], sectionPoints=DEFAULT, 
        variables=('RF2',))
     
    #LET'S CREATE JOBS 
    mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
        explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
        memory=90, memoryUnits=PERCENTAGE, model='Model-%d' %(q), modelPrint=OFF, 
        multiprocessingMode=DEFAULT, name='Job-%d' %(q) , nodalOutputPrecision=SINGLE, 
        numCpus=1, queue=None, scratch='', type=ANALYSIS, userSubroutine='', 
        waitHours=0, waitMinutes=0)
    mdb.jobs['Job-%d' %(q)].writeInput()
    mdb.jobs['Job-%d' %(q) ].submit(consistencyChecking=OFF)    
    mdb.jobs['Job-%d' %(q) ].waitForCompletion()

#Property of Not Real Engineering 
# Author: Shank S. Kulkarni