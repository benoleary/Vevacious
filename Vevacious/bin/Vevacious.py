# This Python program was automatically generated by Vevacious, a program
# released under the GNU General Public License, and, as such, this program
# is also under the GNU General Public License. Please see the
# README.Vevacious.txt file accompanying the Vevacious C++ code for a full
# list of files, brief documentation, and further details on the license.
#
# Vevacious authors: Ben O'Leary (benjamin.oleary@gmail.com)
#                    Florian Staub (florian.staub@googlemail.com)
#                    Jose' Eliel Camargo Molina (elielx@gmail.com)
#                    Werner Porod (porod@physik.uni-wuerzburg.de)
#
#      Copyright 2012 Ben O'Leary, Florian Staub,
#                     Jose' Eliel Camargo Molina, Werner Porod
#
from __future__ import division
import math
import numpy
import numpy.linalg
import minuit
import VevaciousParameterDependent as VPD
import VevaciousTreeLevelExtrema as VTE

warningMessages = []
namesOfVevs = VPD.namesOfVevs
numberOfVevs = len( namesOfVevs )
stepSize = ( 1.0 / VPD.energyScale )
rollingToleranceSquared = ( VPD.rollingTolerance * VPD.rollingTolerance )

effectivePotentialFunction = VPD.LoopCorrectedPotential
# one could replace VPD.LoopCorrectedPotential
# with VPD.TreeLevelPotential for a
# tree-level analysis, for example.

potentialAtVevOrigin = ( VPD.energyScaleFourth
                         * VPD.FunctionFromDictionary( effectivePotentialFunction,
                            VPD.vevOrigin ) )

# MINUIT's hesse() function assumes that it is already at a
# minimum, but we need to check whether it actually stopped at a
# saddle point, so we need to work out the Hessian matrix ourselves.
def NumericalHessian( inputFunction, vevValues ):
    returnHessian = [ [ 0.0 for vevIndexOne in range( numberOfVevs ) ]
                                 for vevIndexTwo in range( numberOfVevs ) ]
    for vevIndexOne in range( numberOfVevs ):
        firstPoint = vevValues.copy()
        firstStepSize = stepSize
        firstPoint[ namesOfVevs[ vevIndexOne ] ] -= firstStepSize
        firstDifference = ( VPD.FunctionFromDictionary( inputFunction,
                                                        vevValues )
                            - VPD.FunctionFromDictionary( inputFunction,
                                                          firstPoint ) )
        for vevIndexTwo in range( vevIndexOne, numberOfVevs ):
            secondPoint = vevValues.copy()
            secondStepSize = stepSize
            secondPoint[ namesOfVevs[ vevIndexTwo ] ] += secondStepSize
            doubleOffset = secondPoint.copy()
            doubleOffset[ namesOfVevs[ vevIndexOne ] ] -= firstStepSize
            returnHessian[ vevIndexOne ][ vevIndexTwo
                         ] = ( ( VPD.FunctionFromDictionary( inputFunction,
                                                             secondPoint )
                               - VPD.FunctionFromDictionary( inputFunction,
                                                             doubleOffset )
                                 - firstDifference )
                               / ( firstStepSize * secondStepSize ) )
            returnHessian[ vevIndexTwo ][ vevIndexOne
                         ] = returnHessian[ vevIndexOne ][ vevIndexTwo ]
    return returnHessian

pointsToTry = VTE.pointsToTry
if ( 0 >= len( pointsToTry ) ):
    warningMessage = "no tree-level extrema found!"
    warningMessages.append( warningMessage )
    print( warningMessage )

foundMinima = []
foundSaddles = []

minuitObject = minuit.Minuit( effectivePotentialFunction )
# There is a possibility of MINUIT trying to roll off to infinity but
# generating an overflow error in calculating the potential before throwing
# an exception. Here limits are set on the VEVs of one thousand times
# the energy scale, which is probably way beyond where the potential should
# be trusted anyway.
for vevVariable in minuitObject.values:
    minuitObject.limits[ vevVariable ] = ( -1.0E+3, 1.0E+3 )


def VevsHaveCorrectSigns( vevValues ):
#    print( "checking signs of " + VPD.UserVevsAsMathematica( vevValues ) )
    returnBool = True
    for positiveVev in VPD.vevsTakenPositive:
        if ( -stepSize > vevValues[ positiveVev ] ):
            returnBool = False
#    print( "returning " + str( returnBool ) )
    return returnBool

def SteepestDescent( vevValues ):
    eigensystemOfHessian = numpy.linalg.eigh( NumericalHessian(
                                                effectivePotentialFunction,
                                                              vevValues ) )
    mostNegativeEigenvalueValue = 0.0
    mostNegativeEigenvalueIndex = 0
    for eigenvalueIndex in range( len( eigensystemOfHessian[ 0 ] ) ):
        if ( mostNegativeEigenvalueValue > eigensystemOfHessian[ 0 ][
                                                       eigenvalueIndex ] ):
            mostNegativeEigenvalueValue = eigensystemOfHessian[ 0 ][
                                                          eigenvalueIndex ]
            mostNegativeEigenvalueIndex = eigenvalueIndex
    return [ mostNegativeEigenvalueValue,
             eigensystemOfHessian[ 1 ][ :, mostNegativeEigenvalueIndex ] ]

def TryToMinimize( vevValues ):
#    print( "trying to minimize" + VPD.UserVevsAsMathematica( vevValues ) )
    global minuitObject
    global foundMinima
    global foundSaddles
    try:
        minuitObject.values = vevValues.copy()
        minuitObject.migrad()
        foundExtremum = { 'potentialDepth': minuitObject.fval,
                          'depthError': abs( minuitObject.edm ),
                          'vevValues': minuitObject.values.copy() }
        candidateDescent = SteepestDescent( minuitObject.values )
        if ( 0.0 > candidateDescent[ 0 ] ):
            foundSaddles.append( [ foundExtremum,
                                   list( candidateDescent[ 1 ] ) ] )
        else:
            foundMinima.append( foundExtremum )
    except minuit.MinuitError as minuitError:
        warningMessage = ( "PyMinuit had problems starting at "
                           + VPD.UserVevsAsMathematica( vevValues )
                           + " [minuit.MinuitError: "
                           + str( minuitError )
                           + "]. PyMinuit stopped at "
                          + VPD.UserVevsAsMathematica( minuitObject.values )
                           + " with relative depth "
               + str( ( VPD.energyScaleFourth
                        * VPD.FunctionFromDictionary( VPD.LoopCorrectedPotential,
                                                    minuitObject.values ) )
                      - potentialAtVevOrigin )
                            + " at one-loop level and "
                 + str( VPD.energyScaleFourth
                        * VPD.FunctionFromDictionary( VPD.TreeLevelPotential,
                                                    minuitObject.values ) )
                            + " at tree level."
               + " Minuit's estimate of how much deeper it should go is "
                 + str( VPD.energyScaleFourth
                        * minuitObject.edm ) )
        warningMessages.append( warningMessage )
        print( warningMessage )

def DisplacePoint( pointDictionary, displacementList, scaleFactor ):
    returnDictionary = pointDictionary.copy()
    for vevIndex in range( numberOfVevs ):
        returnDictionary[ namesOfVevs[ vevIndex ] ] += ( scaleFactor
                                           * displacementList[ vevIndex ] )
    return returnDictionary

for vevValueSet in pointsToTry:
    if ( VevsHaveCorrectSigns( vevValueSet ) ):
        TryToMinimize( vevValueSet )

for saddleSplitNudge in VPD.saddleSplitNudges:
    if ( 0 < len( foundSaddles ) ):
        print( "PyMinuit had to be nudged off "
               + str( len( foundSaddles ) ) + " saddle point(s)." )
        pointsToTry = []
        for saddlePoint in foundSaddles:
            pointsToTry.append( DisplacePoint(
                                         saddlePoint[ 0 ][ 'vevValues' ],
                                                       saddlePoint[ 1 ], 
                                                       saddleSplitNudge ) )
            pointsToTry.append( DisplacePoint(
                                         saddlePoint[ 0 ][ 'vevValues' ],
                                                       saddlePoint[ 1 ], 
                                                      -saddleSplitNudge ) )
        foundSaddles = []
        for vevValueSet in pointsToTry:
            if ( VevsHaveCorrectSigns( vevValueSet ) ):
                TryToMinimize( vevValueSet )

globalMinimumCandidates = foundMinima

# This bit of code adds in any remaining saddle points to
# globalMinimumCandidates as well as setting up warning messages:
if ( 0 != len( foundSaddles ) ):
    warningMessage = ( str( len( foundSaddles ) )
                       + " extremum/a with at least one descending or"
                     + " flat direction remained after all nudging:" )
    for saddlePoint in foundSaddles:
        warningMessage += "\n " + VPD.UserVevsAsMathematica(
                                        saddlePoint[ 0 ][ 'vevValues' ] )
        globalMinimumCandidates.append( saddlePoint[ 0 ] )
    warningMessages.append( warningMessage )
    print( warningMessage )

# If the input vacuum is the global minimum, actionValue is set to -1.0.
# Non-negative values of actionValue indicate the current upper bound on
# the action after the last approximation. It starts stupidly high
# (the current age of the Universe corresponds to an action of about 400)
# so that the first requested bounding estimate will be calculated.
actionNeedsToBeCalculated = False
stabilityVerdict = "error"
actionValue = 1000000.0
actionType = "error"
tunnelingTime = -2.0

if ( 0 >= len( globalMinimumCandidates ) ):
    warningMessage = "no 1-loop-level extrema found!"
    warningMessages.append( warningMessage )
    print( warningMessage )

# The result is assumed long-lived metastable unless found otherwise.
stabilityVerdict = "long-lived"
givenInputAsArray = VPD.VevDictionaryToArray( VPD.inputVevsPoint )
minuitObject.values = VPD.inputVevsPoint.copy()
foundSaddles = []
try:
    minuitObject.migrad()
except minuit.MinuitError as minuitError:
    warningMessage = ( "PyMinuit had problems starting at input VEVs! "
                     + VPDUserVevsAsMathematica( vevValues )
                       + " [minuit.MinuitError: "
                       + str( minuitError )
                       + "]. PyMinuit stopped at "
                       + VPD.UserVevsAsMathematica( minuitObject.values )
                       + " with relative depth "
               + str( ( VPD.energyScaleFourth
                        * VPD.FunctionFromDictionary( VPD.LoopCorrectedPotential,
                                                    minuitObject.values ) )
                      - potentialAtVevOrigin )
                            + " at one-loop level and "
                 + str( VPD.energyScaleFourth
                        * VPD.FunctionFromDictionary( VPD.TreeLevelPotential,
                                                    minuitObject.values ) )
                        + " at tree level."
               + " Minuit's estimate of how much deeper it should go is "
                 + str( VPD.energyScaleFourth
                        * minuitObject.edm ) )
    warningMessages.append( warningMessage )
    print( warningMessage )
if ( 0 != len( foundSaddles ) ):
    warningMessage = "PyMinuit rolled from input VEVs to a saddle point!"
    warningMessages.append( warningMessage )
    print( warningMessage )
rolledInputAsDictionary = minuitObject.values.copy()
# Occasionally degenerate minima with different signs, that are equivalent
# under a gauge transformation, appear, & by the nature of the algorithm,
# MINUIT might roll closer to the sign-flipped minimum than how close it
# rolls to the input minimum. In such cases, the (extremely long) tunneling
# time to the gauge-equivalent minimum would be calculated. To avoid this,
# the input minimum is conservatively taken to have a depth equal to that
# found by MINUIT minus twice MINUIT's error.
rolledInputDepth = minuitObject.fval - abs( 2.0 * minuitObject.edm )
# The input VEV point is taken as the global minimum to begin with:
globalMinimumPointAsDictionary = rolledInputAsDictionary.copy()
globalMinimumDepthError = abs( minuitObject.edm )
globalMinimumDepthValue = rolledInputDepth
rolledInputAsArray = VPD.VevDictionaryToArray( rolledInputAsDictionary )
rollDistanceSquared = numpy.sum( ( rolledInputAsArray
                                   - givenInputAsArray )**2 )
if ( ( rollingToleranceSquared * numpy.sum( givenInputAsArray**2 ) )
     < rollDistanceSquared ):
    warningMessage = (
                  "PyMinuit rolled quite far from the input VEVs! (from "
             + str( VPD.UserVevsAsMathematica( VPD.inputVevsPoint ) ) + " to "
             + str( VPD.UserVevsAsMathematica( minuitObject.values ) ) + ")" )
    warningMessages.append( warningMessage )
    print( warningMessage )
rolledInputLengthSquared = numpy.sum( rolledInputAsArray**2 )

outputFile = open( "./Vevacious_loop-corrected_results.txt", "w" )

for globalMinimumCandidate in globalMinimumCandidates:
    outputFile.write( str( globalMinimumCandidate[ 'potentialDepth' ]
                           * VPD.energyScaleFourth ) + ", "
                      + VPD.UserVevsAsMathematica(
                      globalMinimumCandidate[ 'vevValues' ] ) + "\n" )
    if ( globalMinimumDepthValue
         > globalMinimumCandidate[ 'potentialDepth' ] ):
        actionNeedsToBeCalculated = True
        globalMinimumPointAsDictionary = globalMinimumCandidate[
                                                     'vevValues' ].copy()
        globalMinimumDepthValue = globalMinimumCandidate[
                                                       'potentialDepth' ]
        globalMinimumDepthError = globalMinimumCandidate[ 'depthError' ]

outputFile.close()

numericallyDegenerateCandidates = []
for globalMinimumCandidate in globalMinimumCandidates:
    if ( ( globalMinimumCandidate[ 'depthError' ]
           + globalMinimumDepthError )
         > ( globalMinimumCandidate[ 'potentialDepth' ]
             - globalMinimumDepthValue ) ):
       numericallyDegenerateCandidates.append(
                                            globalMinimumCandidate.copy() )

shortestDistanceSquared = -1.0
for closeEnoughMinimum in numericallyDegenerateCandidates:
    if ( rolledInputDepth > closeEnoughMinimum[ 'potentialDepth' ] ):
        closeEnoughMinimumAsArray = VPD.VevDictionaryToArray(
                                      closeEnoughMinimum[ 'vevValues' ] )
        distanceSquaredToInput = numpy.sum( ( closeEnoughMinimumAsArray
                                              - rolledInputAsArray )**2 )
# shortestDistanceSquared is negative before the 1st
# candidate minimum has been checked for its distance from the input
# VEVs, so this 1st "if" ensures that the 1st candidate minimum
# is taken as the best minimum so far automatically.
        if ( shortestDistanceSquared < 0.0 ):
            shortestDistanceSquared = ( distanceSquaredToInput + 1.0 )
        if ( distanceSquaredToInput < shortestDistanceSquared ):
            shortestDistanceSquared = distanceSquaredToInput
            globalMinimumPointAsDictionary = closeEnoughMinimum[
                                                     'vevValues' ].copy()

globalMinimumPointAsArray = VPD.VevDictionaryToArray(
                                           globalMinimumPointAsDictionary )

# if the global minimum is sufficiently far away from the input VEVs that
# it is unlikely that is was just that MINUIT did not roll exactly to where
# there should have been a minimum at the input VEVs...
distanceSquaredFromInputToGlobalMinimum = numpy.sum( (
                      globalMinimumPointAsArray - rolledInputAsArray )**2 )

# If the input vacuum is the global minimum, actionValue is set to -1.0.
# Non-negative values of actionValue indicate the current upper bound on
# the action after the last approximation. It starts stupidly high
# (the current age of the Universe corresponds to an action of about 400)
# so that the first requested bounding estimate will be calculated.
if ( ( rollingToleranceSquared * rolledInputLengthSquared )
     > distanceSquaredFromInputToGlobalMinimum ):
    stabilityVerdict = "stable"
    actionValue = -1.0
    tunnelingTime = -1.0
    actionType = "unnecessary"
    actionNeedsToBeCalculated = False

# The resolution of the tunneling path needs to be set
# (low-ish by default for speed):
tunnelingResolution = 20

if ( actionNeedsToBeCalculated
     and ( ( 0.0 < VPD.directActionBound )
           or ( 0.0 < VPD.deformedActionBound ) ) ):
    firstStepPoint = ( ( rolledInputAsArray
                         * ( 1.0 - ( 1.0 / tunnelingResolution ) ) )
                       + ( globalMinimumPointAsArray
                           * ( 1.0 / tunnelingResolution ) ) )
    firstStepDepth = VPD.FunctionFromArray( effectivePotentialFunction,
                                                firstStepPoint )
    if ( rolledInputDepth >= firstStepDepth ):
        actionType = "barrier_smaller_than_resolution"
        actionValue = 0.0
        tunnelingTime = 0.0
        stabilityVerdict = "short-lived"
        actionNeedsToBeCalculated = False
        warningMessage = ( "Energy barrier from input VEVs to global"
                        + " minimum thinner than resolution of tunneling"
                           + " path!" )
        warningMessages.append( warningMessage )
        print( warningMessage )

if ( actionNeedsToBeCalculated
     and ( ( 0.0 < VPD.directActionBound )
           or ( 0.0 < VPD.deformedActionBound ) ) ):
    import sys
    sys.path.append( VPD.pathToCosmotransitions )
    import pathDeformation as CPD

    arrayOfArrays = numpy.array( [ globalMinimumPointAsArray.copy(),
                                   rolledInputAsArray.copy() ] )

    def PotentialFromArray( pointAsArray ):
        return VPD.FunctionFromArray( effectivePotentialFunction,
                                                             pointAsArray )

    def PotentialFromMatrix( arrayOfArrays ):
        if ( ( numberOfVevs, ) == arrayOfArrays.shape ):
            return PotentialFromArray( arrayOfArrays )
        elif ( ( len( arrayOfArrays ), numberOfVevs )
               == arrayOfArrays.shape ):
            returnArray = numpy.zeros( len( arrayOfArrays ) )
            for whichIndex in range( len( arrayOfArrays ) ):
                returnArray[ whichIndex ] = PotentialFromArray(
                                              arrayOfArrays[ whichIndex ] )
            return returnArray
        else:
            return None

    def GradientFromArray( pointAsArray ):
        potentialAtPoint = PotentialFromArray( pointAsArray )
        gradientArray = numpy.zeros( len( pointAsArray ) )
        for whichField in range( len( pointAsArray ) ):
            displacedPoint = pointAsArray.copy()
            displacedPoint[ whichField ] += stepSize
            gradientArray[ whichField ] = (
                                     ( PotentialFromArray( displacedPoint )
                                       - potentialAtPoint ) / stepSize )
        return gradientArray

    def GradientFromMatrix( arrayOfArrays ):
        if ( ( numberOfVevs, ) == arrayOfArrays.shape ):
            return GradientFromArray( arrayOfArrays )
        elif ( ( len( arrayOfArrays ), numberOfVevs )
               == arrayOfArrays.shape ):
            returnMatrix = arrayOfArrays.copy()
            for whichIndex in range( len( arrayOfArrays ) ):
                returnMatrix[ whichIndex ] = GradientFromArray(
                                              arrayOfArrays[ whichIndex ] )
            return returnMatrix
        else:
            return None

    if ( ( 0.0 < VPD.directActionBound )
         and ( actionValue > VPD.directActionBound ) ):
        quickTunneler = CPD.fullTunneling( V = PotentialFromMatrix,
                                           dV = GradientFromMatrix,
                                           phi = arrayOfArrays,
                                           quickTunneling = False,
                                           npoints = tunnelingResolution )
        quickTunneler.tunnel1D( xtol = 1e-4, phitol = 1e-6 )
        actionValue = quickTunneler.findAction()
        actionType = "direct_path_bound"
        if( actionValue < VPD.directActionBound ):
            stabilityVerdict = "short-lived"
            actionNeedsToBeCalculated = False

    if ( actionNeedsToBeCalculated
         and ( 0.0 < VPD.deformedActionBound )
         and ( actionValue > VPD.deformedActionBound ) ):
        fullTunneler = CPD.fullTunneling( V = PotentialFromMatrix,
                                          dV = GradientFromMatrix,
                                          phi = arrayOfArrays,
                                          quickTunneling = False,
                                          npoints = tunnelingResolution )
# setting a maximum of 20 path deformation iterations before giving up on
# finding the optimal path may seem defeatist, but I have rarely seen it
# converge if it hasn't within the first few iterations. an action is still
# calculated, though it may not be the minimum action possible.
        fullTunneler.run( maxiter = 20 )
        actionValue = fullTunneler.findAction()
        actionType = "full_deformed_path"
        if ( actionValue < VPD.deformedActionBound ):
            stabilityVerdict = "short-lived"
            actionNeedsToBeCalculated = False

# No matter if there were serious errors or not, an output file is written:
outputFile = open( VPD.outputFile, "w" )
outputFile.write( "<Vevacious_result>\n"
                  + "  <reference version=\"0.3.0\" citation=\"[still unpublished]\" />\n"
             + "  <stability> " + stabilityVerdict + " </stability>\n"
                  + "  <global_minimum   relative_depth=\""
                      + str( ( globalMinimumDepthValue * VPD.energyScaleFourth ) - potentialAtVevOrigin ) + "\" "
                      + VPD.UserVevsAsXml( globalMinimumPointAsDictionary )
                      + " />\n  <input_minimum   relative_depth=\""
                      + str( ( rolledInputDepth * VPD.energyScaleFourth ) - potentialAtVevOrigin ) + "\" "
                      + VPD.UserVevsAsXml( rolledInputAsDictionary )
                      + " />" )
if ( 0.0 <= actionValue ):
# We don't want an overflow error when calculating the tunneling time.
    if ( 1000.0 < actionValue ):
        actionValue = 1000.0
    tunnelingTime = ( math.exp( 0.25 * actionValue ) * 1.0e-43 )
    if ( 1000000.0 < tunnelingTime ):
        tunnelingTime = 1000000.0
# The tunneling time is capped at a million times the current age of the
# known Universe.
# This code assumes that the age of the Universe is (10^(44))/TeV, and that
# the solitonic solution factor (Sidney Coleman's "A") is (0.1 TeV)^4.
outputFile.write( "\n  <lifetime  action_calculation=\"" + actionType
                  + "\" > " + str( tunnelingTime ) + " </lifetime>" )
# Each warning is printed as an XML element:
for warningMessage in warningMessages:
    outputFile.write( "\n  <warning>\n  " + warningMessage
                      + "\n  </warning>" )
outputFile.write( "\n</Vevacious_result>\n" )
outputFile.close()
