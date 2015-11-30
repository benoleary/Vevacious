/*
 * PotentialMinimizer.cpp
 *
 *  Created on: Oct 5, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of Vevacious, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.Vevacious.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#include "../include/PotentialMinimizer.hpp"

namespace Vevacious
{
  std::string const PotentialMinimizer::overallFactorAttributeName( "factor" );
  std::string const
  PotentialMinimizer::subtractionConstantAttributeName( "constant" );
  std::string const
  PotentialMinimizer::spinTypeAttributeName( "spin" );
  std::string const
  PotentialMinimizer::defaultSubtractionConstantString( "1.5" );


  PotentialMinimizer::PotentialMinimizer() :
    sarahParser( false ),
    attributeFinder(),
    pythonCode(),
    polynomialPartOfPotential( "" ),
    factoredMassSquaredMatrices(),
    massesSquaredFunctions(),
    stringParser(),
    argumentsWithoutBrackets( "" ),
    argumentsWithBrackets( "()" )
  {
    pythonCode.precision( 17 );
    stringParser.precision( 17 );
    // HOM4PS2 operates at 17 digits of precision.
    // SARAH should prepare all the mass-squared matrices (even for fermions)
    // as XML elements. there is no set number of mass matrices.
  }

  PotentialMinimizer::~PotentialMinimizer()
  {
    sarahParser.closeFile();
  }


  void PotentialMinimizer::addMassSquaredMatrix(
                      std::map< std::string, std::string > const& attributeMap,
                                         std::string const& massSquaredMatrix )
  {
    attributeFinder
    = attributeMap.find( subtractionConstantAttributeName );
    std::string subtractionConstant( defaultSubtractionConstantString );
    if( attributeMap.end() != attributeFinder )
    {
      subtractionConstant.assign( BOL::StringParser::trimFromFrontAndBack(
                                                       attributeFinder->second,
                                                                    "\"\'" ) );
    }
    attributeFinder = attributeMap.find( spinTypeAttributeName );
    std::string spinType( "" );
    if( attributeMap.end() != attributeFinder )
    {
      spinType.assign( BOL::StringParser::trimFromFrontAndBack(
                                                       attributeFinder->second,
                                                                    "\"\'" ) );
    }
    attributeFinder = attributeMap.find( overallFactorAttributeName );
    if( attributeMap.end() != attributeFinder )
    {
      factoredMassSquaredMatrices.push_back( StringQuadruple() );
      factoredMassSquaredMatrices.back().first.first.assign(
              BOL::StringParser::trimFromFrontAndBack( attributeFinder->second,
                                                       "\"\'" ) );
      factoredMassSquaredMatrices.back().first.second.assign(
                                                         subtractionConstant );
      factoredMassSquaredMatrices.back().second.first.assign( spinType );
      factoredMassSquaredMatrices.back().second.second.assign(
                                                           massSquaredMatrix );
      SarahInterpreter::removeNewlinesFrom(
                            factoredMassSquaredMatrices.back().second.second );
    }
  }

  void PotentialMinimizer::prepareLoopCorrections(
                                           SarahInterpreter& sarahInterpreter )
  {
    massesSquaredFunctions.clear();
    std::string massSquaredMatrix( "" );
    std::string overallFactor( "" );
    std::string subtractionConstant( "" );
    std::string massSquaredFunctionName( "" );
    BOL::VectorlikeArray< std::string > massTerms;
    for( std::vector< StringQuadruple >::iterator
         whichQuadruple( factoredMassSquaredMatrices.begin() );
         factoredMassSquaredMatrices.end() > whichQuadruple;
         ++whichQuadruple )
    {
      massTerms.clearEntries();
      BOL::StringParser::parseByChar( sarahInterpreter(
                                               whichQuadruple->second.second ),
                                      massTerms,
                                      ',' );
      stringParser.clear();
      stringParser.str( "" );
      stringParser
      << "MassesSquared" << ( massesSquaredFunctions.size() + 1 );
      massSquaredFunctionName.assign( stringParser.str() );
      pythonCode <<
"# The masses-squared from this function are in GeV^2:\n"
"def " << massSquaredFunctionName << argumentsWithBrackets << ":\n"
"    massSquaredElementArray = numpy.array( [ ";
      for( int whichTerm( 0 );
           massTerms.getSize() > whichTerm;
           ++whichTerm )
      {
        if( 0 < whichTerm )
        {
          pythonCode << ", ";
        }
        pythonCode << " ( "
        << BOL::StringParser::trimFromFrontAndBack( massTerms[ whichTerm ],
                                 BOL::StringParser::whitespaceAndNewlineChars )
        << " )";
      }
      pythonCode << " ] )\n"
"    massSquaredMatrix = massSquaredElementArray.reshape( math.sqrt(\n"
"                                       massSquaredElementArray.size ), -1 )\n"
"    return numpy.linalg.eigvalsh( massSquaredMatrix )\n"
"\n";
      stringParser
      << ", " << whichQuadruple->first.first
      << ", " << whichQuadruple->first.second
      << ", \"" << whichQuadruple->second.first << "\" ]";
      massesSquaredFunctions.push_back( "[ " + stringParser.str() );
      // all the corrections have to be noted for later.
    }
    pythonCode <<
"\n"
"\n"
"MassesSquareds = [ ";
    for( std::vector< std::string >::iterator
         whichMassesSquared( massesSquaredFunctions.begin() );
         massesSquaredFunctions.end() > whichMassesSquared;
         ++whichMassesSquared )
    {
      if( whichMassesSquared > massesSquaredFunctions.begin() )
      {
        pythonCode << ",\n                   ";
      }
      pythonCode << *whichMassesSquared;
    }
    pythonCode << " ]\n"
"\n"
"def LoopCorrectedPotential( " << argumentsWithoutBrackets
<<                                              ", temperatureValue = 0.0 ):\n"
"    loopCorrections = 0.0\n"
"    for MassesSquared in MassesSquareds:\n"
"        loopCorrections += MassSquaredCorrections( MassesSquared[ 0 ]"
<<                                               argumentsWithBrackets << ",\n"
"                                                   MassesSquared[ 1 ],\n"
"                                                   MassesSquared[ 2 ] )\n"
"    return ( PolynomimalPartOfPotential" << argumentsWithBrackets << "\n"
"             + ( loopFactor * loopCorrections ) )\n"
"\n"
"\n"
"def LoopAndThermalCorrectedPotential( " << argumentsWithoutBrackets
<<                                                    ", temperatureValue ):\n"
"    loopCorrections = 0.0\n"
"    thermalCorrections = 0.0\n"
"    temperatureSquared = ( temperatureValue * temperatureValue )\n"
"    for MassesSquared in MassesSquareds:\n"
"        massesSquared = MassesSquared[ 0 ]" << argumentsWithBrackets << "\n"
"        loopCorrections += MassSquaredCorrections( massesSquared,\n"
"                                                   MassesSquared[ 1 ],\n"
"                                                   MassesSquared[ 2 ] )\n"
"        if ( temperatureValue > 0.0 ):\n"
"            adjustedOverallFactor = MassesSquared[ 1 ]\n"
"            if ( MassesSquared[ 3 ] is \"vector\" ):\n"
"                adjustedOverallFactor = ( ( 2.0 * adjustedOverallFactor )\n"
"                                          / 3.0 )\n"
"            thermalCorrections += ThermalCorrections( massesSquared,\n"
"                                                     adjustedOverallFactor,\n"
"                                                      temperatureSquared )\n"
"    return ( PolynomimalPartOfPotential" << argumentsWithBrackets << "\n"
"             + ( loopFactor * loopCorrections )\n"
"             + ( thermalFactor * temperatureSquared * temperatureSquared\n"
"                 * thermalCorrections ) )\n"
"\n"
"\n"
"def AbsLoopAndThermalCorrectedPotential( " << argumentsWithoutBrackets
<<                                                    ", temperatureValue ):\n"
"    loopCorrections = 0.0\n"
"    thermalCorrections = 0.0\n"
"    temperatureSquared = ( temperatureValue * temperatureValue )\n"
"    for MassesSquared in MassesSquareds:\n"
"        massesSquared = MassesSquared[ 0 ]" << argumentsWithBrackets << "\n"
"        loopCorrections += MassSquaredCorrections( massesSquared,\n"
"                                                   MassesSquared[ 1 ],\n"
"                                                   MassesSquared[ 2 ] )\n"
"        if ( temperatureValue > 0.0 ):\n"
"            adjustedOverallFactor = MassesSquared[ 1 ]\n"
"            if ( MassesSquared[ 3 ] is \"vector\" ):\n"
"                adjustedOverallFactor = ( ( 2.0 * adjustedOverallFactor )\n"
"                                          / 3.0 )\n"
"            thermalCorrections += AbsThermalCorrections( massesSquared,\n"
"                                                     adjustedOverallFactor,\n"
"                                                       temperatureSquared )\n"
"    return ( PolynomimalPartOfPotential" << argumentsWithBrackets << "\n"
"             + ( loopFactor * loopCorrections )\n"
"             + ( thermalFactor * temperatureSquared * temperatureSquared\n"
"                 * thermalCorrections ) )\n"
"\n"
"\n"
"def FloorLoopAndThermalCorrectedPotential( " << argumentsWithoutBrackets
<<                                                    ", temperatureValue ):\n"
"    loopCorrections = 0.0\n"
"    thermalCorrections = 0.0\n"
"    temperatureSquared = ( temperatureValue * temperatureValue )\n"
"    for MassesSquared in MassesSquareds:\n"
"        massesSquared = MassesSquared[ 0 ]" << argumentsWithBrackets << "\n"
"        loopCorrections += MassSquaredCorrections( massesSquared,\n"
"                                                   MassesSquared[ 1 ],\n"
"                                                   MassesSquared[ 2 ] )\n"
"        if ( temperatureValue > 0.0 ):\n"
"            adjustedOverallFactor = MassesSquared[ 1 ]\n"
"            if ( MassesSquared[ 3 ] is \"vector\" ):\n"
"                adjustedOverallFactor = ( ( 2.0 * adjustedOverallFactor )\n"
"                                          / 3.0 )\n"
"            thermalCorrections += FloorThermalCorrections( massesSquared,\n"
"                                                     adjustedOverallFactor,\n"
"                                                       temperatureSquared )\n"
"    return ( PolynomimalPartOfPotential" << argumentsWithBrackets << "\n"
"             + ( loopFactor * loopCorrections )\n"
"             + ( thermalFactor * temperatureSquared * temperatureSquared\n"
"                 * thermalCorrections ) )\n"
"\n"
"\n"
"# PotentialScaler class:\n"
"\n"
"class PotentialScaler:\n"
"    \"\"\"\n"
"    This class exists solely to provide a function that is a scaled version\n"
"    of the function given to the constructor for MINUIT, and a set of\n"
"    functions using arrays that the CosmoTransitions objects want to use.\n"
"    Well, it also changes the function beyond a hypersurface of radiuszn"
"    fieldLimit: within the hypersurface, the function is just scaled, but\n"
"    outside, the function is taken to be the value it had on the surface at\n"
"    the same angular co-ordinates, plus the square of the difference of the\n"
"    square of the Euclidean length of the field configuration from the\n"
"    square of fieldLimit. I.e., if the field configuration vector is longer\n"
"    than fieldLimit, it is scaled to be of length fieldLimit, and then the\n"
"    difference is its original length squared minus fieldLimit squared, and\n"
"    this difference is itself squared and added to the value of the\n"
"    potential evaluated at the scaled field configuration vector.\n"
"    \"\"\"\n"
"\n"
"    def __init__( self,\n"
"                  FunctionToScale,\n"
"                  fieldLimit,\n"
"                  temperatureValue = 0.0 ):\n"
"        self.FunctionToScale = FunctionToScale\n"
"        self.functionAtOrigin = FunctionToScale( ";
    std::vector< char > const&
    vevsOrderedBySolution( sarahInterpreter.getVevsOrderedBySolutions() );
    for( unsigned int whichVev( 0 );
        vevsOrderedBySolution.size() > whichVev;
        ++whichVev )
    {
      pythonCode << "0.0, ";
    }
    pythonCode << "0.0 )\n"
"        self.fieldLimitSquared = ( inverseScale * fieldLimit )**2\n"
"        self.temperatureValue = temperatureValue\n"
"\n"
"\n"
"    def SetMaximumConfigurationLength( self, fieldLimit ):\n"
"        self.fieldLimitSquared = ( inverseScale * fieldLimit )**2\n"
"\n"
"\n"
"    def ScaledFunctionFromScaledArguments( self,\n"
"                                      " << argumentsWithoutBrackets << " ):\n"
"        configurationLengthSquared = ( ";
    for( std::vector< char >::const_iterator
         whichName( vevsOrderedBySolution.begin() );
         vevsOrderedBySolution.end() > whichName;
         ++whichName )
    {
      if( vevsOrderedBySolution.begin() != whichName )
      {
        pythonCode << " + ";
      }
      pythonCode << "( " << *whichName << " )**2";
    }
    pythonCode << " )\n"
"        lengthSquaredBeyondCap = 0.0\n"
"        if ( configurationLengthSquared > self.fieldLimitSquared ):\n"
"            lengthSquaredBeyondCap = ( configurationLengthSquared\n"
"                                       - self.fieldLimitSquared )\n"
"            scalingRatio = math.sqrt( self.fieldLimitSquared\n"
"                                      / configurationLengthSquared )\n";
    for( std::vector< char >::const_iterator
         whichName( vevsOrderedBySolution.begin() );
         vevsOrderedBySolution.end() > whichName;
         ++whichName )
    {
      pythonCode << "            " << *whichName << " = ( scalingRatio * "
      << *whichName << " )\n";
    }
    pythonCode << "        functionValue = ( self.FunctionToScale( ";
    for( std::vector< char >::const_iterator
         whichName( vevsOrderedBySolution.begin() );
         vevsOrderedBySolution.end() > whichName;
         ++whichName )
    {
      pythonCode
      << *whichName << " = ( " << *whichName
      << " * energyScale ),\n                                                ";
    }
    pythonCode << "temperatureValue = self.temperatureValue )\n"
"                          - self.functionAtOrigin\n"
"                          + lengthSquaredBeyondCap**2 )\n"
"        return ( inverseScaleFourthed\n"
"                 * functionValue )\n"
"\n"
"    def PotentialFromArray( self, pointAsArray ):\n"
"        return FunctionFromArray( self.FunctionToScale,\n"
"                                  pointAsArray,\n"
"                                  self.temperatureValue )\n"
"\n"
"    def PotentialFromMatrix( self, arrayOfArrays ):\n"
"        if ( ( numberOfFields, ) == arrayOfArrays.shape ):\n"
"            return self.PotentialFromArray( arrayOfArrays )\n"
"        elif ( ( len( arrayOfArrays ), numberOfFields )\n"
"               == arrayOfArrays.shape ):\n"
"            returnArray = numpy.zeros( len( arrayOfArrays ) )\n"
"            for whichIndex in range( len( arrayOfArrays ) ):\n"
"                returnArray[ whichIndex ] = self.PotentialFromArray(\n"
"                                              arrayOfArrays[ whichIndex ] )\n"
"            return returnArray\n"
"        else:\n"
"            return None\n"
"\n"
"    def GradientFromArray( self, pointAsArray ):\n"
"        potentialAtPoint = self.PotentialFromArray( pointAsArray )\n"
"        gradientArray = numpy.zeros( len( pointAsArray ) )\n"
"        for whichField in range( len( pointAsArray ) ):\n"
"            displacedPoint = pointAsArray.copy()\n"
"            displacedPoint[ whichField ] += numericalStepSize\n"
"            gradientArray[ whichField ] = ( ( self.PotentialFromArray(\n"
"                                                           displacedPoint )\n"
"                                              - potentialAtPoint )\n"
"                                            / numericalStepSize )\n"
"        return gradientArray\n"
"\n"
"    def GradientFromMatrix( self, arrayOfArrays ):\n"
"        if ( ( numberOfFields, ) == arrayOfArrays.shape ):\n"
"            return self.GradientFromArray( arrayOfArrays )\n"
"        elif ( ( len( arrayOfArrays ), numberOfFields )\n"
"               == arrayOfArrays.shape ):\n"
"            returnMatrix = arrayOfArrays.copy()\n"
"            for whichIndex in range( len( arrayOfArrays ) ):\n"
"                returnMatrix[ whichIndex ] = self.GradientFromArray(\n"
"                                              arrayOfArrays[ whichIndex ] )\n"
"            return returnMatrix\n"
"        else:\n"
"            return None\n"
"\n"
"# End of PotentialScaler class\n"
"\n";
  }

} /* namespace Vevacious */
