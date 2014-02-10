/*
 * SarahInterpreter.cpp
 *
 *  Created on: Feb 7, 2013
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2013 Ben O'Leary
 *
 *      This file is part of Vevacious, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.Vevacious.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#include "../include/SarahInterpreter.hpp"

namespace Vevacious
{
  std::string const SarahInterpreter::imaginaryUnitChars( "iIjJ" );

  SarahInterpreter::SarahInterpreter() :
    slhaConverter( NULL ),
    treeBlockPrefix( "TREE" ),
    loopBlockPrefix( "LOOP" ),
    vevRenamer(),
    argumentString( "" ),
    returnString( "" ),
    nextCharPosition( 0 ),
    lastCharPosition( 0 )
  {
    // just an initialization list.
  }

  SarahInterpreter::~SarahInterpreter()
  {
    delete slhaConverter;
  }


  void SarahInterpreter::setPositiveVevs( std::string const& positiveVevs )
  {
    positiveInternalVevs.assign( "[ " );
    BOL::AsciiXmlParser xmlParser;
    if( xmlParser.loadString( positiveVevs ) )
    {
      while( xmlParser.readNextElement() )
      {
        if( xmlParser.currentElementNameMatches( "taken_positive" ) )
        {
          returnString.assign( xmlParser.getCurrentElementContent() );
          BOL::StringParser::substituteCharacterWith( returnString,
                                                      ',',
                                                      ' ' );
          BOL::VectorlikeArray< std::string > positiveUserVevs;
          BOL::StringParser::parseByChar( returnString,
                                          positiveUserVevs,
                                BOL::StringParser::whitespaceAndNewlineChars );
          for( int whichVev( 0 );
               positiveUserVevs.getSize() > whichVev;
               ++whichVev )
          {
            if( 0 < whichVev )
            {
              positiveInternalVevs.append( ", " );
            }
            positiveInternalVevs.append( "\'" );
            positiveInternalVevs.append( 1,
               vevRenamer.getInternalVevName( positiveUserVevs[ whichVev ] ) );
            positiveInternalVevs.append( "\'" );
          }
          break;
        }
      }
    }
    positiveInternalVevs.append( " ]" );
  }

  std::string const&
  SarahInterpreter::operator()( std::string const& inputString,
                                bool const useTreeRatherThanLoop,
                                bool const writingTadpole )
  // this returns a reference to a string where all SLHA block elements in
  // inputString have been replaced by their values and all the user-defined
  // VEV names in inputString have been replaced by their internal
  // single-char names. if emptyIfHasSlhaZero is true and any of the SLHA
  // substitutions was for zero, a reference to an empty string is returned
  // instead. if emptyIfHasImaginaryUnit is true and any chars from
  // charsDenotingImaginaryUnit were in the string after accounting for SLHA
  // substitutions, a reference to an empty string is returned instead.
  // useTreeRatherThanLoop determines which prefix is prepended to SLHA block
  // names as a 1st replacement attempt, before the block name without a prefix
  // is used if the prepended block name is not found for the appropriate
  // index.
  {
    char currentChar( '?' );
    if( useTreeRatherThanLoop )
    {
      argumentString.assign( treeBlockPrefix );
    }
    else
    {
      argumentString.assign( loopBlockPrefix );
    }
    std::string const& slhaReplacedString( (*slhaConverter)( inputString,
                                                             argumentString,
                                                             writingTadpole,
                                                             "0.0" ) );
    std::string slhaAndVevReplacedString( vevRenamer.replaceUserVevNames(
                                                        slhaReplacedString ) );
    if( writingTadpole
        &&
        ( std::string::npos
          != slhaAndVevReplacedString.find_first_of( imaginaryUnitChars ) ) )
    {
      return VevRenamer::emptyString;
    }
    nextCharPosition = slhaAndVevReplacedString.find_first_of( "^IiEe" );
    if( std::string::npos == nextCharPosition )
    {
      returnString.assign( slhaAndVevReplacedString );
      return returnString;
    }
    lastCharPosition = 0;
    returnString.assign( slhaAndVevReplacedString,
                         0,
                         nextCharPosition );
    while( std::string::npos != nextCharPosition )
    {
      currentChar = slhaAndVevReplacedString[ nextCharPosition ];
      if( '^' == currentChar )
      {
        if( writingTadpole )
        {
          returnString.append( "^" );
        }
        else
        {
          returnString.append( "**" );
        }
      }
      else if( ( 'i' == currentChar )
               ||
               ( 'I' == currentChar ) )
      {
        returnString.append( "(1.0j)" );
      }
      else
      {
        if( writingTadpole
            &&
            ( slhaAndVevReplacedString.size() > ( nextCharPosition + 1 ) )
            &&
            !( ( '+' == slhaAndVevReplacedString[ nextCharPosition + 1 ] )
               ||
               ( '-' == slhaAndVevReplacedString[ nextCharPosition + 1 ] ) ) )
          // if 'e' or 'E' is not immediately followed by '+' or '-'...
        {
          // then we assume that the exponent should be a positive number.
          returnString.append( "E+" );
        }
        else
        {
          returnString.append( 1,
                               currentChar );
        }
      }
      lastCharPosition = ( nextCharPosition + 1 );
      nextCharPosition = slhaAndVevReplacedString.find_first_of( "^IiEe",
                                                            lastCharPosition );
      returnString.append( slhaAndVevReplacedString,
                           lastCharPosition,
                           ( nextCharPosition - lastCharPosition ) );
    }
    return returnString;
  }

} /* namespace Vevacious */
