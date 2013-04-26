/*
 * VevaciousRunner.hpp
 *
 *  Created on: Mar 5, 2013
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2013 Ben O'Leary
 *
 *      This file is part of Vevacious, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.Vevacious.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef VEVACIOUSRUNNER_HPP_
#define VEVACIOUSRUNNER_HPP_

#include <fstream>
#include <stdexcept>
#include <sstream>
#include "BOLlib/include/ArgumentParser.hpp"
#include "BOLlib/include/AsciiXmlParser.hpp"
#include "BOLlib/include/UsefulStuff.hpp"
#include "BOLlib/include/WaitingOnSubprocessExecutor.hpp"
#include "SarahInterpreter.hpp"
#include "TadpoleSolver.hpp"
#include "PotentialMinimizer.hpp"

namespace Vevacious
{
  // this class wraps all the Vevacious functionality together.
  class VevaciousRunner
  {
  public:
    VevaciousRunner( std::string const& modelFilename,
                     std::string const& hom4ps2Directory,
                     std::string const& homotopyType );
    VevaciousRunner( BOL::ArgumentParser& argumentParser );
    ~VevaciousRunner();

    void setResultsFilename( std::string const& resultsFilename )
    { this->resultsFilename.assign( resultsFilename ); }
    void setImaginaryPartTolerance( double imaginaryTolerance )
    { this->imaginaryTolerance = imaginaryTolerance; }
    void setImaginaryPartTolerance( std::string const& imaginaryTolerance )
    { this->imaginaryTolerance
      = BOL::StringParser::stringToDouble( imaginaryTolerance ); }
    void setMaximumMinuitNudgesOffSaddlePoints( int maximumSaddleSplits )
    { this->maximumSaddleSplits = maximumSaddleSplits; }
    void setMaximumMinuitNudgesOffSaddlePoints(
                                       std::string const& maximumSaddleSplits )
    { this->maximumSaddleSplits
      = BOL::StringParser::stringToInt( maximumSaddleSplits ); }
    void setMinuitRollingTolerance( double rollingTolerance )
    { this->rollingTolerance = rollingTolerance; }
    void setMinuitRollingTolerance( std::string const& rollingTolerance )
    { this->rollingTolerance
      = BOL::StringParser::stringToDouble( rollingTolerance ); }
    void setLifetimeForQuarticGuess( double quarticLifetimeBound )
    { quarticActionBound = actionFromLifetime( quarticLifetimeBound ); }
    void setLifetimeForQuarticGuess( std::string const& quarticLifetimeBound )
    { setLifetimeForQuarticGuess(
                 BOL::StringParser::stringToDouble( quarticLifetimeBound ) ); }
    void setPathToCosmotransitions( std::string const& pathToCosmotransitions )
    { this->pathToCosmotransitions.assign( pathToCosmotransitions ); }
    void setLifetimeForDirectPath( double directLifetimeBound )
    { directActionBound = actionFromLifetime( directLifetimeBound ); }
    void setLifetimeForDirectPath( std::string const& directLifetimeBound )
    { setLifetimeForDirectPath(
                  BOL::StringParser::stringToDouble( directLifetimeBound ) ); }
    void setLifetimeForDeformedPath( double deformedLifetimeBound )
    { deformedActionBound = actionFromLifetime( deformedLifetimeBound ); }
    void setLifetimeForDeformedPath( std::string const& deformedLifetimeBound )
    { setLifetimeForDeformedPath(
                BOL::StringParser::stringToDouble( deformedLifetimeBound ) ); }
    void appendTreeLevelExtrema( std::string const& slhaFilename,
                                 std::string const solutionsFilename
                                 = "VevaciousTreeLevelExtrema.py" );
    void overwriteTreeLevelExtrema( std::string const& slhaFilename,
                                    std::string const solutionsFilename
                                    = "VevaciousTreeLevelExtrema.py" );
    void writeCalculatedTreeLevelExtrema( std::string const& outputFilename );
    void prepareParameterDependentPython( std::string const& slhaFilename,
                                          std::string const pythonFilename
                                          = "VevaciousParameterDependent.py" );
    void runPython( std::string const pythonFilename = "Vevacious.py" );
    // this checks to see if there is a file called pythonFilename, then, if it
    // exists, runs it with the terminal's "python". if it doesn't exist, then
    // the default file with name given by defaultPythonFilename is written &
    // run.
    bool runPython( std::string const pythonFilename,
                    int const deformationPatience );
    // this does the same as runPython( std::string const pythonFilename ), but
    // killing the Python subprocess after deformationPatience seconds if it
    // takes that long. it returns true if the subprocess finished without
    // having to be killed, false if the subprocess took too long.
    void appendResultsToSlha( std::string const& slhaFilename,
                   bool const properFormatRatherThanSspReadable = true ) const;
    // this takes the XML-format results from the file called resultsFilename,
    // & appends the results & warnings in custom SLHA blocks to the end of the
    // file name slhaFilename.


  protected:
    static double const lifetimeFactor;
    static std::string const pointsToTry;
    static std::string const resultsFilenameVariableName;
    static std::string const maximumSaddleSplitsVariableName;
    static std::string const rollingToleranceVariableName;
    static std::string const quarticActionBoundVariableName;
    static std::string const pathToCosmotransitionsVariableName;
    static std::string const directActionBoundVariableName;
    static std::string const deformedActionBoundVariableName;
    static std::string const defaultPythonFilename;
    static std::string const treeLevelExtremaHeader;
    static BOL::StringParser const slhaIndexMaker;
    static BOL::StringParser const slhaDoubleMaker;

    static std::string slhaDoubleFromQuotedString(
                                             std::string const& quotedString );

    SarahInterpreter sarahInterpreter;
    TadpoleSolver tadpoleSolver;
    PotentialMinimizer potentialMinimizer;
    std::stringstream treeLevelExtrema;
    bool firstWriteOfExtrema;
    std::string resultsFilename;
    double imaginaryTolerance;
    int maximumSaddleSplits;
    double rollingTolerance;
    double quarticActionBound;
    std::string pathToCosmotransitions;
    double directActionBound;
    double deformedActionBound;

    void calculateTreeLevelExtrema( std::string const& slhaFilename );
    void writeDefaultPythonProgram() const;
    std::string
    prepareToRunPythonCommand( std::string const& pythonFilename ) const;
    double actionFromLifetime( double lifetimeValue ) const;
  };





  inline void
  VevaciousRunner::appendTreeLevelExtrema( std::string const& slhaFilename,
                                          std::string const solutionsFilename )
  {
    BOL::UsefulStuff::runSystemCommand( "rm " + solutionsFilename + "c" );
    calculateTreeLevelExtrema( slhaFilename );
    std::ofstream outputFile( solutionsFilename.c_str(),
                              std::ios::app );
    if( !(outputFile.good()) )
    {
      throw std::runtime_error( "could not open \"" + solutionsFilename
                                + "\" to write tree-level extrema." );
    }
    if( firstWriteOfExtrema )
    {
      outputFile << treeLevelExtremaHeader;
      firstWriteOfExtrema = false;
    }
    outputFile << treeLevelExtrema.str();
    outputFile.close();
  }

  inline void
  VevaciousRunner::overwriteTreeLevelExtrema( std::string const& slhaFilename,
                                          std::string const solutionsFilename )
  {
    BOL::UsefulStuff::runSystemCommand( "rm " + solutionsFilename + "c" );
    calculateTreeLevelExtrema( slhaFilename );
    std::ofstream outputFile( solutionsFilename.c_str() );
    if( !(outputFile.good()) )
    {
      throw std::runtime_error( "could not open \"" + solutionsFilename
                                + "\" to write tree-level extrema." );
    }
    outputFile << treeLevelExtremaHeader << treeLevelExtrema.str();
    outputFile.close();
  }

  inline void VevaciousRunner::runPython( std::string const pythonFilename )
  // this checks to see if there is a file called pythonFilename, then, if it
  // exists, runs it with the terminal's "python". if it doesn't exist, then
  // the default file with name given by defaultPythonFilename is written &
  // run.
  {
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "VevaciousRunner::runPython( \"" << pythonFilename << "\" ) called.";
    std::cout << std::endl;/**/

    BOL::UsefulStuff::runSystemCommand(
                                 prepareToRunPythonCommand( pythonFilename ) );
  }

  inline bool VevaciousRunner::runPython( std::string const pythonFilename,
                                    int const deformationPatience )
   // this does the same as runPython( std::string const pythonFilename ), but
   // killing the Python subprocess after deformationPatience seconds if it
   // takes that long. it returns true if the subprocess finished without
   // having to be killed, false if the subprocess took too long.
   {
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "VevaciousRunner::runPython( \"" << pythonFilename << "\", "
    << deformationPatience << " ) called.";
    std::cout << std::endl;/**/

     std::string systemCommand( prepareToRunPythonCommand( pythonFilename ) );
     if( ( 0 < deformationPatience )
         &&
         ( 0.0 < deformedActionBound ) )
     {
       std::string forkedCommand( "bash " + systemCommand );
       BOL::WaitingOnSubprocessExecutor subprocessExecutor( forkedCommand,
                                                            true,
                                               ( 1000 * deformationPatience ) );
       // BOL::WaitingOnSubprocessExecutor takes a time in milliseconds, while
       // deformedPathPatience is in seconds.
       return subprocessExecutor.forkAndExecvAndWait();
     }
     BOL::UsefulStuff::runSystemCommand( systemCommand );
     return true;
   }

  inline std::string VevaciousRunner::slhaDoubleFromQuotedString(
                                              std::string const& quotedString )
  {
    return slhaDoubleMaker.doubleToString( BOL::StringParser::stringToDouble(
                         BOL::StringParser::trimFromFrontAndBack( quotedString,
                                                                  "\'\"" ) ) );
  }

  inline std::string VevaciousRunner::prepareToRunPythonCommand(
                                      std::string const& pythonFilename ) const
  {
    std::string systemCommand( "rm " + pythonFilename + "c" );
    BOL::UsefulStuff::runSystemCommand( systemCommand );
    systemCommand.assign( "python " );
    if( !(BOL::UsefulStuff::fileExists( pythonFilename )) )
    {
      writeDefaultPythonProgram();
      systemCommand.append( defaultPythonFilename );
    }
    else
    {
      systemCommand.append( pythonFilename );
    }
    return systemCommand;
  }

  inline double
  VevaciousRunner::actionFromLifetime( double lifetimeValue ) const
  {
    if( !( 0.0 < lifetimeValue ) )
    {
      return -1.0;
    }
    return ( ( 4.0 * log( lifetimeValue ) ) + lifetimeFactor );
  }

} /* namespace Vevacious */
#endif /* VEVACIOUSRUNNER_HPP_ */
