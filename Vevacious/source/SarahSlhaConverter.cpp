/*
 * SarahSlhaConverter.cpp
 *
 *  Created on: Oct 2, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of Vevacious, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.Vevacious.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#include "../include/SarahSlhaConverter.hpp"

namespace Vevacious
{
  std::string const SarahSlhaConverter::blockKeyword( "SLHA::" );
  char const SarahSlhaConverter::blockIndexOpener( '[' );
  char const SarahSlhaConverter::blockIndexCloser( ']' );
  std::string const SarahSlhaConverter::blockIndexChars( "[ .,0123456789]" );

  SarahSlhaConverter::SarahSlhaConverter( std::string const& slhaFilename ) :
    slhaValues( slhaFilename ),
    blockStart( 0 ),
    blockEnd( 0 ),
    searchString( "" ),
    slhaString( "" ),
    returnString( "" )
  {
    // just an initialization list.
  }

  SarahSlhaConverter::~SarahSlhaConverter()
  {
    // does nothing.
  }


  std::string const&
  SarahSlhaConverter::operator()( std::string const& inputString,
                                  std::string const& prefixString,
                                  bool const zeroesReturnEmptyString,
                                  std::string const emptyElementReplacement )
  // this sets returnString to be the same as inputString with all SLHA block
  // elements in inputString replaced by their values, then returns a
  // reference to returnString. however, if no SLHA block elements were
  // found, a reference to stringToConvert is returned instead. if
  // emptyIfHasSlhaZero is true and any of the SLHA substitutions was for
  // zero, a reference to emptyString is returned instead.
  {
    blockStart = inputString.find( blockKeyword );
    // blockStart now gives the position in inputString of the 1st char of the
    // 1st instance of blockKeyword in inputString.
    if( std::string::npos == blockStart )
      // if there were no blocks to replace, the original string is returned.
    {
      return inputString;
    }
    // otherwise, set up returnString & once it is ready, return a reference to
    // it.
    returnString.assign( inputString.substr( 0,
                                             blockStart ) );
    while( std::string::npos != blockStart )
    {
      // returnString is now all of inputString up to the 1st instance of
      // blockKeyword.
      blockStart += blockKeyword.size();
      // blockStart is now the position in inputString of the 1st char after
      // the 1st instance of blockKeyword in inputString.
      blockEnd = inputString.find( blockIndexCloser,
                                   inputString.find( blockIndexOpener,
                                                     blockStart ) );
      // blockEnd is now the position of the last char of the SLHA block
      // keyword to replace.
      if( std::string::npos == blockEnd )
        // if the string had a malformed SLHA block keyword, a reference to
        // emptyString is returned.
      {
        return VevRenamer::emptyString;
      }
      searchString.assign( prefixString + inputString.substr( blockStart,
                                                 ( blockEnd - blockStart ) ) );
      // the string to look up in slhaValues is the substring immediately
      // following the last found instance of blockKeyword, up until the last
      // char of the chars giving the indices of the block element, with
      // prefixString prepended for the 1st check, & without prefixString
      // prepended for the check that follows if slhaValues did not find the
      // prepended string.
      slhaString.assign( slhaValues.withMap( searchString ) );
      if( slhaString.empty() )
      {
        slhaString.assign( slhaValues.withMap( inputString.substr( blockStart,
                                               ( blockEnd - blockStart ) ) ) );
      }
      if( zeroesReturnEmptyString
          &&
          ( slhaString.empty()
            ||
            ( 0.0 == BOL::StringParser::stringToDouble( slhaString ) ) ) )
      {
        // if the block element was 0 or not mentioned in the SLHA file, the
        // term is considered to be multiplied by zero & thus stringToConvert
        // is set to be an empty string if zeroesReturnEmptyString was set to
        // true.
        return VevRenamer::emptyString;
      }
      // otherwise, the string from the SLHA file is substituted in:
      returnString.append( "(" );
      if( slhaString.empty() )
      {
        returnString.append( emptyElementReplacement );
      }
      else
      {
        returnString.append( slhaString );
      }
      returnString.append( ")" );
      // now we look for further SLHA block keywords:
      blockStart = inputString.find( blockKeyword,
                                     (++blockEnd) );
      // blockEnd is incremented because it was at the position of the last
      // character of the block string, rather than at the 1st character after.
      returnString.append( inputString,
                           blockEnd,
                           ( blockStart - blockEnd ) );
    }
    return returnString;
  }

} /* namespace Vevacious */
