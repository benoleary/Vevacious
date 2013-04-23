/*
 * TadpoleSolver.hpp
 *
 *  Created on: Oct 2, 2012
 *      Authors: José Eliel Camargo (elielx@gmail.com),
 *               Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 José Eliel Camargo, Ben O'Leary
 *
 *      This file is part of Vevacious, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.Vevacious.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef TADPOLESOLVER_HPP_
#define TADPOLESOLVER_HPP_

#include <set>
#include <vector>
#include <fstream>
#include <cstdio>
#include <sstream>
#include <string>
#include <stdexcept>
#include <limits.h>
#include "BOLlib/include/ArgumentParser.hpp"
#include "BOLlib/include/UsefulStuff.hpp"
#include "SLHA.hpp"
#include "SarahInterpreter.hpp"
#include "VevRenamer.hpp"

namespace Vevacious
{
  // this class writes tadpole equations from SARAH with SLHA block values
  // substituted in as an input file for HOM4PS2, then calls HOM4PS2, then
  // parses the output of HOM4PS2.
  class TadpoleSolver
  {
  public:
    TadpoleSolver( std::string const& hom4ps2Directory,
                   std::string const& homotopyType );
    ~TadpoleSolver();

    void setEquations( std::string const& tadpolesString );
    void solveTadpoles( SarahInterpreter& sarahInterpreter,
                        double const zeroImaginaryTolerance = 0.000001 );
    // this writes an input file for HOM4PS2 based on the tadpole equations
    // provided by SARAH with the values from slhaValues, runs HOM4PS2, and
    // parses its output.
    void writeTreeLevelSolutions( std::string const& outputFilename );
    // this writes the solutions from HOM4PS2 in a Mathematica-friendly format
    // in outputFilename.
    std::set< std::vector< long double > > const&
    getPurelyRealSolutionSets() const;


  protected:
    static std::string const relativeHom4ps2InputFilename;
    static std::string const relativeHom4ps2OutputFilename;

    std::string const hom4ps2Directory;
    std::string const homotopyType;
    std::string const absoluteHom4ps2InputFilename;
    std::string const absoluteHom4ps2OutputFilename;

    std::stringstream parsingStream;
    BOL::VectorlikeArray< std::string > tadpoleEquations;
    std::string currentConvertedTadpole;
    std::vector< std::string > convertedTadpoles;
    BOL::CommentedTextParser tadpoleSolutionsFile;

    std::vector< std::pair< long double, long double > > complexSolutions;
    std::vector< char > variableNames;
    std::set< std::vector< long double > > purelyRealSolutionSets;
    int purelyRealSolutionSetCount;
    std::set< std::vector< std::pair< long double, long double > > >
    notPurelyRealSolutionSets;
    int notPurelyRealSolutionSetCount;

    void writeSpecificTadpoles( SarahInterpreter& sarahInterpreter );
    // this creates an input file for HOM4PS2 using tadpoleEquations &
    // sarahSlhaConverter.
    void runHomotopy();
    // this makes a note of the current working directory in
    // originalWorkingDirectory, then changes directory to hom4ps2Directory,
    // then runs HOM4PS2 with absoluteHom4ps2InputFilename as input, with
    // homotopyType as well, then returns to originalWorkingDirectory.
    void parseSolutions( SarahInterpreter& sarahInterpreter,
                         double const zeroImaginaryTolerance );
    // this parses the output of HOM4PS2 into 2 set of sets of solutions of the
    // system of tadpole equations: 1 set where each set has only real
    // solutions, & the other set has the rest of the sets. the solution sets
    // are stored as long doubles or pairs of long doubles.
    std::string longDoubleForMathematica( long double longDouble );
    // this returns a string that is "((" + (longDouble with a precision of 17
    // digits) + "))", replacing 'E' or 'e' with ")*10^(" (hence the "((" &
    // "))", so it is consistent).
    void convertTadpoleEquation( std::string const& tadpoleString,
                                 SarahInterpreter& sarahInterpreter );
    // this puts tadpoleString with its SLHA values substituted in into
    // currentConvertedTadpole.
  };




  inline void TadpoleSolver::setEquations( std::string const& tadpolesString )
  {
    tadpoleEquations.clearEntries();
    std::string editedTadpolesString( tadpolesString + " " );
    SarahInterpreter::removeNewlinesFrom( editedTadpolesString );
    BOL::StringParser::parseByChar( editedTadpolesString,
                                    tadpoleEquations,
                                    ';' );
    // if there are n tadpole equations, each ending with ';', tadpoleEquations
    // is now size ( n + 1 ). the added " " is to ensure that the last entry
    // of tadpoleEquations is an invalid equation, so can be ignored. (it
    // breaks otherwise if any whitespace would get in after the ';' of the
    // last tadpole.)
    if( 1 > tadpoleEquations.getSize() )
    {
      std::cout << std::endl
      << "Error! could not find any tadpole equations!";
      std::cout << std::endl;
      throw std::runtime_error( "no tadpole equations" );
    }
  }

  inline void TadpoleSolver::solveTadpoles(
                                           SarahInterpreter& sarahInterpreter,
                                          double const zeroImaginaryTolerance )
  // this writes an input file for HOM4PS2 based on the tadpole equations
  // provided by SARAH with the values from slhaValues, runs HOM4PS2, and
  // parses its output.
  {
    // the function names should explain the process well enough...
    writeSpecificTadpoles( sarahInterpreter );
    runHomotopy();
    parseSolutions( sarahInterpreter,
                    zeroImaginaryTolerance );
  }

  inline std::set< std::vector< long double > > const&
  TadpoleSolver::getPurelyRealSolutionSets() const
  {
    return purelyRealSolutionSets;
  }

  inline std::string TadpoleSolver::longDoubleForMathematica(
                                                       long double longDouble )
  // this returns a string that is "((" + (longDouble with a precision of 17
  // digits) + "))", replacing 'E' or 'e' with ")*10^(" (hence the "((" &
  // "))", so it is consistent).
  {
    std::string returnString( "" );
    parsingStream.clear();
    parsingStream.str( "" );
    parsingStream << "((" << longDouble;
    parsingStream >> returnString;
    size_t positionOfExponentChar( returnString.find_first_of( "Ee" ) );
    if( std::string::npos != positionOfExponentChar )
    {
      returnString[ positionOfExponentChar ] = ')';
      returnString.insert( positionOfExponentChar,
                           "*10^(" );
    }
    returnString.append( "))");
    return returnString;
  }

} /* namespace Vevacious */
#endif /* TADPOLESOLVER_HPP_ */
