/*
 * SarahSlhaConverter.hpp
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

#ifndef SARAHSLHACONVERTER_HPP_
#define SARAHSLHACONVERTER_HPP_

#include <string>
#include <set>
#include <stdexcept>
#include "SLHA.hpp"
#include "VevRenamer.hpp"

namespace Vevacious
{
  // this class replaces substrings corresponding to an SLHA block element as
  // written by SARAH in strings.
  class SarahSlhaConverter
  {
  public:
    SarahSlhaConverter( std::string const& slhaFilename );
    ~SarahSlhaConverter();


    std::string const& operator()( std::string const& inputString,
                     std::string const& prefixString = VevRenamer::emptyString,
                                   bool const zeroesReturnEmptyString = false,
                           std::string const emptyElementReplacement = "0.0" );
    // this sets returnString to be the same as inputString with all SLHA block
    // elements in inputString replaced by their values, then returns a
    // reference to returnString. however, if no SLHA block elements were
    // found, a reference to stringToConvert is returned instead. if
    // emptyIfHasSlhaZero is true and any of the SLHA substitutions was for
    // zero, a reference to emptyString is returned instead.
    double getScale() const;


  protected:
    static std::string const blockKeyword;
    static char const blockIndexOpener;
    static char const blockIndexCloser;
    static std::string const blockIndexChars;

    LHPC::SlhaSimplisticInterpreter slhaValues;
    size_t blockStart;
    size_t blockEnd;
    std::string searchString;
    std::string slhaString;
    std::string returnString;
    double gaugeScale;
  };




  inline double SarahSlhaConverter::getScale() const
  {
    return gaugeScale;
  }

} /* namespace Vevacious */

#endif /* SARAHSLHACONVERTER_HPP_ */
