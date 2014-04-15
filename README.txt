/*****************************************************************************\
 * README.txt                                                                *
 *                                                                           *
 *  Created on: Oct 8, 2012                                                  *
 *      Authors: Ben O'Leary (benjamin.oleary@gmail.com)                     *
 *      Copyright 2012 Ben O'Leary                                           *
 *                                                                           *
 *      This file is part of Vevacious, a program designed to find the       *
 *      configuration of vacuum expectation values (VEVs) for the global     *
 *      minimum of a quantum field theory potential energy function.         *
 *                                                                           *
 *      Vevacious is free software: you can redistribute it and/or modify it *
 *      under the terms of the GNU General Public License as published by    *
 *      the Free Software Foundation, either version 3 of the License, or    *
 *      (at your option) any later version.                                  *
 *                                                                           *
 *      Vevacious is distributed in the hope that it will be useful, but     *
 *      WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU     *
 *      General Public License for more details.                             *
 *                                                                           *
 *      You should have received a copy of the GNU General Public License    *
 *      (in HoTTMiLC/GNU_public_license.txt ) along with Vevacious. If not,  *
 *      see <http://www.gnu.org/licenses/>.                                  *
 *      A full list of the files of Vevacious is at the end of this file.    *
\*****************************************************************************/

 Now that the legalese preamble is out of the way, the description of the code
 and how to use it can begin!

 Vevacious is a program written in C++ to try to find the global minimum of a
 loop-corrected potential. It is designed to be used in conjuction with SARAH
 by Florian Staub and utilizes HOM4PS2 by Tsung-Lin Lee, Tien-Yien Li, and
 Chih-Hsiung Tsai, PyMinuit by Jim Pivarski, implementing the MINUIT
 algorithms by Fred James, and CosmoTransitions, by Carroll (Max) Wainwright.
 
 Vevacious is an implementation of the process of attempting to find the global
 minimum of a potential as decided upon by the collaboration of
 José Eliel Camargo (elielx@gmail.com),
 Ben O'Leary (benjamin.oleary@gmail.com),
 Werner Porod (porod@physik.uni-wuerzburg.de), and
 Florian Staub (florian.staub@googlemail.com). This README file was written by
 Ben O'Leary.
 
 A manual describing the physics and process of Vevacious is available from
 http://arxiv.org/abs/1307.1477 and there is a quickstart guide available from
 http://vevacious.hepforge.org/ following the appropriate link. An rough
 version of the installation instructions is also contained in this README
 file, below the description of the process.
 
 Vevacious runs in 2 stages: first it prepares the tree-level tadpole equation
 system and then solves it with HOM4PS2, and second, it uses the tree-level
 solutions as starting points for minimizing the 1-loop-level potential with
 PyMinuit, including, if a deeper minimum than that desired was found,
 calculating a tunneling time or an upper bound using CosmoTransitions.
 
 The first part writes a system of tadpole equations for a given set of
 parameters (provided in SLHA format, presumably from a version of SPheno
 created by SARAH) for a given model (created by SARAH) using the general form
 of the tree-level tadpole equations in a file provided by SARAH, in the format
 required by HOM4PS2. it then runs HOM4PS2 and parses its output. At this
 stage, the tree-level extrema can be printed out as a file in a
 Mathematica-friendly format (which I hope is also reasonably human-readable).
 In principle the homotopy continutation algorithm finds _all_ the extrema of
 the tree-level potential, but finite-precision and finite-step-size issues
 mean that any algorithm will sometimes miss some extrema. The user should be
 aware of this.

 The second part takes the results of HOM4PS2 as the starting points for
 minimizing the loop-corrected potential with PyMinuit. It does this by writing
 a program in Python with the tree-level potential and all the loop corrections
 for the specific values of the Lagrangian parameters found in the SLHA file
 from the general forms taken from a file provided by SARAH for the given
 model. This potential is then minimized by PyMinuit starting from each unique
 extremum found by HOM4PS2 (ignoring solutions with complex roots, because
 SARAH already writes any complex fields as 2 real degrees of freedom, so only
 real roots are valid solutions of the provided system of tadpole equations).
 The deepest minimum thus found is written in the output file in XML, which,
 again, I hope is reasonably human-readable. (The other minimization results
 are available by default in Vevacious_loop-corrected_results.txt in a
 Mathematica-friendly format.)
 There is no guarantee that the loop corrections do _not_ introduce new minima
 beyond the number of minima that the tree-level potential has. By default,
 Vevacious will move the MINUIT object off saddle points in the steepest
 direction and its mirror, so will find minima that develop when a tree-level
 minimum becomes a one-loop-level saddle point (such as the origin for SPS1a'
 in a certain renormalization scheme).
 There is also no guarantee that the MINUIT minimization starting from the
 tree-level extrema will find as many minima as the tree-level potential has,
 regardless of whether the loop corrections introduce more minima. The user
 should be aware that this program tries to find the global minimum of the
 loop-corrected potential, but it could fail to find it if the loop corrections
 are sufficiently large.

 So, if the above limitations are not a problem, the following instructions
 should allow a successful installation:
 1) Download SARAH. The files are available at
    http://sarah.hepforge.org/
    (link last checked 2013-08-29). The installation is as simple as unpacking
    the gzipped tarball.
 2) Download HOM4PS2. The files are available at
    http://www.math.nsysu.edu.tw/~leetsung/works/HOM4PS_soft_files/HOM4PS_Linux.htm
    (link last checked 2013-08-29). The installation is as simple as unpacking
    the gzipped tarball.
 3) Ensure that Python is installed. PyMinuit requires at least version 2.4 or
    later. I shouldn't have to get into any specifics of how to install Python
    here... Internet search engines are your friends.
 4) Download and install PyMinuit. The instructions are available at
    http://code.google.com/p/pyminuit/wiki/HowToInstall
    (link last checked 2013-08-29). The installation involves downloading and
    compiling one of the C++ implementations of MINUIT (such as the standalone
    version 1.7.9 or the CERN ROOT version, both available from
    http://lcgapp.cern.ch/project/cls/work-packages/mathlibs/minuit/
    (link last checked 2013-08-29)), then running the PyMinuit setup script
    (giving the path where the C++ MINUIT code was _built_, *not* installed -
    it needs the .o object files rather than the .a library file...). The
    LD_LIBRARY_PATH and PYTHONPATH environment variables then need to be set.
 5) Download CosmoTransitions. The files are available at
    http://chasm.uchicago.edu/cosmotransitions/
    (link last checked 2013-08-29). The installation is as simple as unpacking
    the zipped file.
    (CosmoTransitions was previously at http://chasm.ucsc.edu/cosmotransitions/
    but has since moved. Hosting on HepForge at some point in the future has
    been suggested, but as of 2013-08-29, it is not there.)
    WARNING! CosmoTransitions v1.0.2 no longer works with the most recent
    versions of SciPy and NumPy! Replacing "integrate.inf" with "numpy.inf"
    in the CosmoTransitions code fixes this.
 6) Download and compile the LesHouchesParserClasses (LHPC) C++ library. The
    files are available at
    http://www.hepforge.org/downloads/lhpc
    (link last checked 2013-08-29) or
    https://github.com/benoleary/LesHouchesParserClasses_CPP
    (link last checked 2013-08-29). The installation should just be
    downloading, unzipping, and then running make.
 7) Compile Vevacious. The Makefile should be edited to have the correct paths
    to the header files and library file of LHPC (these are
    /path/to/unzipped/LHPC/include/ and /path/to/unzipped/LHPC/lib/ by
    default). LHPC version 0.8.4 or higher is required.

 Once Vevacious is installed, it can be run as is, as long as the paths are set
 correctly. Various parameters for each execution can be specified: for the
 list, see the Vevacious/bin/VevaciousInitialization.xml file, where the
 various tags <XYZ> can be over-ridden by being given as arguments
 (e.g. --XYZ=ABC will over-ride whatever value is given between <XYZ> and
 </XYZ> in the .ini file). A different initialization file can be given with
 the option --input=other_filename as well (and direct arguments will over-ride
 any given in that file too). The XML file shows what paths must be specified.

 An executable to perform batch runs on a set of input SLHA files with the same
 model file is provided as VevaciousBatch.exe and takes the same initialization
 file, though it ignores any --slha_file/<slha_file> option, instead looking
 for <input_dir> and <output_dir> (as usual, they can be over-ridden by
 command-line arguments, such as --input_dir=/path/to/input/dir/ and so on).
 Each file in input_dir will be taken as an input SLHA file, and will be
 copied to output_dir with VEVACIOUSRESULTS appended to it once Vevacious has
 finished with that point. While Vevacious is working on the parameter point,
 a placeholder file with the extension ".placeholder" is left in output_dir.
 This allows multiple cores (each with their OWN COPY of HOM4PS2 and its ENTIRE
 DIRECTORY - this is really important as HOM4PS2 uses temporary files for
 intermediate results, so 2 Vevacious processes running 2 HOM4PS2 processes is
 likely to cause problems as one HOM4PS2 process uses the temporary file
 generated by the other process) to work through input_dir, without
 over-writing the work of other cores. For example, core1 and core2 look at
 example/path/input/ and find example1.spc and example2.spc. If core1 gets
 there 1st, it makes example/path/output/example1.spc.placeholder, and when it
 is finished, it copies example1.spc to example/path/output/ and appends the
 VEVACIOUSRESULTS block, and deletes the placeholder. Meanwhile, core2 looks at
 example1, but finds that example/path/output/example1.spc.placeholder already
 exists, so it moves on to example2. Later, another core comes along, and for
 example finds that example/path/output/example1.spc already exists and
 example/path/output/example2.spc.placeholder also already exists, so it does
 not try to run example1.spc or example2.spc.
 
 The batch executable will probably crash if any other files or any
 subdirectories are in input_dir. Sorry.


CHANGELOG:
 * 15th April 2014: version 1.1.00beta12
 - Fixed very minor bug in default Python code that was incorrectly comparing
   massSquaredMagnitude to ( 1.0E-6 * inverseScaleSquared ) instead of
   ( 1.0E-6 / inverseScaleSquared ). It was not a problem in the sense that it
   was a check to avoid division by zero and worked even though it was wrong.
   It was however possibly calculating m^4 ln( m^2 / Q^2 ) for very small m
   which would have been a negligible contribution to the corrections. Fixed
   the comparison to just whether massSquaredMagnitude > 1.0 or not.

 * 3rd April 2014: version 1.1.00beta11
 - Fixed bug in default Python code that was incorrectly exluding points based
   on thermal fluctuations based on the dividing an action by the last
   temperature in a list rather than on the temperature at which the action was
   calculated.
 - Fit of direct-path 3-dimensional action in default Python code improved:
   the fit is only made to the action at temperatures above the evaporation
   temperature of the DSB minimum.

 * 25th March 2014: version 1.1.00beta10
 - Added setThermalSurvivalThreshold functions to VevaciousRunner class.
 - Shifted a lot of code from default Vevacious.py to
   VevaciousParameterDependent.py, mainly to do with finding the optimal
   tunneling temperature.
 - Python code to calculate the survival probability against thermal tunneling
   has been improved: some silly bugs leading to misreporting of the result
   have been fixed, and now the result "long-lived_but_thermally_excluded" is
   reserved for exclusion by a cautious lower bound on the integrated decay
   width, while failing to be excluded by this but still being excluded by the
   old, aggressive algorithm is now "long-lived_but_maybe_thermally_excluded".
   The number code reported in the SLHA block VEVACIOUSRESULTS is -2 for
   excluded by the cautious calculation, and -3 for exclusion only by the
   aggressive calculation while not being excluded by the cautious calculation.
 - Default Vevacious.py now also decides the number code for the (0, 0) entry
   of the SLHA block VEVACIOUSRESULTS, which is communicated back to
   Vevacious.exe by an XML attribute.

 * 14th March 2014: version 1.1.00beta9
 - Added flags to turn off thermal calculations and path deformation more
   easily:
   --should_tunnel=True
   --tunnel_thermally=True
   --deform_tunnel_paths=True
   or with "True" replaced by "False" - note that the case is important! It
   must be valid Python, so with uppercase initial letters, and the rest
   lowercase! All values are assumed to be True by default.
 - Included VevaciousBatchRunner.cpp which the Makefile compiles to
   VevaciousBatch.exe as an executable to run Vevacious on all files in a
   directory.
 - Removed erroneous whitespace character appearing in <reference> element
   attribute "version" of output XML produced by default Python.
 - Changed default main program to use the VevaciousRunner constructor which
   takes just a reference to a BOL::ArgumentParser instance.

 * 11th March 2014: version 1.1.00beta8
 - Default Python now generates AbsLoopAndThermalCorrectedPotential and
   FloorLoopAndThermalCorrectedPotential which use the absolute value or 0.0
   for negative masses-squared respectively in the thermal corrections.

 * 4th March 2014: version 1.1.00beta7
 - Default Python now correctly stops if it finds that the point is thermally
   excluded by a direct path when exploring the temperature dependence of the
   direct-path thermal action.
 - Default Python now also considers tunneling from the false vacuum to the
   true vacuum at a given temperature impossible if the true vacuum is now less
   than 0.1 times its length at zero temperature, as a workaround for PyMinuit
   sometimes not rolling the DSB minimum close enough to the field origin when
   it should, breaking the unpatched algorithm for deciding if tunneling is
   still possible.

 * 28th February 2014: version 1.1.00beta6
 - Default parameter-dependent Python now correctly does not deform the path if
   not necessary.

 * 28th February 2014: version 1.1.00beta5
 - Default parameter-dependent Python now correctly reports exclusions in
   output file (was correctly reporting to terminal, but results file was
   wrong for some cases).

 * 17th February 2014: version 1.1.00beta4
 - Default parameter-dependent Python now loads warning message list in
   Vevacious object before trying to minimize the DSB minimum, to stop crashes
   caused by trying to log errors from the minimization of the DSB minimum.

 * 14th February 2014: version 1.1.00beta3
 - Fixed default Vevacious.py to properly exclude field configurations
   that would be the DSB minimum but for numerical effects from the list of
   possible panic vacua.
 - Put function back in to limit MINUIT to a hypercube in field space in
   default Vevacious.py code.
 - Improved fit function for guess at temperature dependence of thermal action
   by direct path in default Python code (now a quadratic divided by the square
   of the difference from the best guess at the critical temperature).

 * 13th February 2014: version 1.1.00beta2
 - Fixed default Vevacious.py to properly report thermal exclusion if found
   during direct path thermal action calculation.
 - Updated bin/Vevacious.py example to the default created by 1.1.00beta2.

 * 11th February 2014: version 1.1.00beta1
 ~ Major update!
 - Added functionality to calculate survival probability against tunneling to
   panic vacua at non-zero temperatures.
   -- Default Python program now
      >> takes tree-level extrema from parsed results of HOM4PS2 as starting
         points for PyMinuit using the zero-temperature loop-corrected
         potential
      >> sorts minima and from those minima deeper than the DSB minimum it
         chooses the nearest as the panic vacuum
      >> checks for exclusion based on tunneling time at zero temperature with
         a direct path between the minima
      >> finds a temperature between 2^-1/2 & 1 times the critical temperature
         at which the panic vacuum becomes less deep than the DSB vacuum (or
         where the DSB vacuum rolls to at the temperature)
      >> checks for exclusion based on the survival probability against thermal
         tunneling along a direct path between the minima at this temperature
         and at half this temperature
      >> fits a guess at the temperature dependence of the direct-path thermal
         action based on the actions at the above temperatures and finds the
         optimal tunneling temperature according to the fitted function
      >> checks for exclusion based on the survival probability against thermal
         tunneling along a direct path between the minima at the estimated
         optimal temperature
      >> checks for exclusion based on tunneling time at zero temperature with
         an optimal deformed path between the minima
      >> checks for exclusion based on the survival probability against thermal
         tunneling along an optimal deformed path between the minima at the
         estimated optimal temperature
      >> stops the above calculation if at any stage the parameter point is
         excluded
 - Most boilerplate Python code is now in VevaciousParameterDependent.py so
   that the main Python program (defaulting to Vevacious.py) can be neater.
 - Model files now need that the mass-squared matrices for vectors have the XML
   element spin="vector" so that thermal corrections can be calculated
   correctly. All bundled example .vin files have this included, and Florian
   will incorporate this into the next update of SARAH.
 - The functions VevaciousRunner::setLifetimeForDirectPath and
   VevaciousRunner::setLifetimeForDeformedPath (both variants of each) have
   been removed in favor of a single function
   VevaciousRunner::setLifetimeThreshold which sets a single threshold (because
   I cannot remember what I was thinking when I set it up so that direct and
   deformed paths would have separate thresholds).
 - Many static strings in VevaciousRunner and PotentialMinimizer have been
   removed in favor of just having their values written into the function that
   writes the Python code.
 - VevaciousRunner::prepareParameterDependentPython and
   VevaciousRunner::writeDefaultPythonProgram have had bits of the boilerplate
   Python code swapped around between them, along with
   PotentialMinimizer::prepareLoopCorrections.
 - Fields are no longer scaled to the energy scale except within the Python, so
   VevRenamer and SarahInterpreter were changed to reflect that.
 ~ Not ready to be 1.1.00 yet, as VevaciousRunner::appendResultsToSlha needs to
   be updated to account for thermal results.

 * 9th October 2013: version 1.0.11
 - Fixed that default Vevacious.py was using the number of spatial dimensions
   for a finite-temperature tunneling time calculation rather than the correct
   zero-temperature calculation number of dimensions.
 - Fixed bug that max_saddle_nudges was being set to 2 if not explicitly given,
   when it should be however not be used to resize saddle_nudges unless
   explicitly given.
 - Changed default homotopy_type from 1 (polyhedral) to 2 (linear) as it seems
   to be faster for SUSY potentials, and thus I reckon that it's better for
   generic QFT potentials.
 - Removed taken_positive from provided model files, as occasionally it leads
   to global minima being missed because HOM4PS2 did not find all the sign
   combinations.

 * 13th September 2013: version 1.0.10
 - Fixed bug when trying to use a relative path for the hom4ps2_dir input.
 - Fixed bug that sometimes a point that has a negative value for a VEV that
   should be positive (as declared by the <taken_positive> element of the
   <input_vevs> element of the model file) was taken as the global minimum at
   tree level, leading to misleading warnings about apparent change from
   metastable to stable going from tree to one loop.

 * 11th September 2013: version 1.0.9
 - Default Python program fixed to correctly find the tree-level global minimum
   for the purposes of checking to see if the basin of attraction of the
   one-loop minima has moved significantly. (The code was not correctly
   indented, leading to stuff being evaluated after a loop rather than during
   the loop.)
 - Added example model file for just Higgs VEVs and stop VEVs (without allowing
   stau VEVs), though just the SARAH-SPhenoMSSM style of SLHA is expected.

 * 9th September 2013: version 1.0.8
 - Example model files where the stop VEVs are allowed to be non-zero have been
   corrected (unfortunately the D-term from SU(3)_c had been generated wrongly
   and this carried through into the mass matrices).

 * 29th August 2013: version 1.0.7
 - The A factor for calculating the tunneling time has changed to be the fourth
   power of renormalization scale as given by the SLHA file rather than the old
   hard-coded (100 GeV)^4 by default, and the default Vevacious.py now
   explicitly calculates the actions from the tunneling time thresholds. This
   now allows the user to edit Vevacious.py to change the calculation of the
   A factor if so desired.
 - CosmoTransitions is now at http://chasm.uchicago.edu/cosmotransitions/ as
   Dr Wainwright has kindly let me know.

 * 21st August 2013: version 1.0.6
 - Makefile fixed so that it works properly (the libraries were in the wrong
   official order, but some compilers don't mind, such as that which was used
   to test 1.0.5).
 - Example model files renamed to clarify (a bit) which require non-standard
   SLHA blocks (by the prefix "SARAH-SPheno"), and which can be used with
   standard flavor-violation-accommodating SLHA2 blocks alone.
 - Example model files that should work with SLHA files adhering just to the
   SLHA1 conventions have been added.
 - Warning! The website at chasm.uscs.edu does not seem to exist any longer, so 
   it appears to be impossible to get CosmoTransitions officially! Please email
   Ben O'Leary if Dr C. Wainwright (the author of CosmoTransitions) cannot help
   you.

 * 1st August 2013: version 1.0.5
 - Makefile now makes ./lib/libVevacious.a as a static library as well, which
   should make it easier to make custom C++ programs that use the Vevacious
   classes.
 - Default Python now also checks to see if a tree-level analysis of metastable
   turned into a 1-loop stable result, & if so, doubles the VEVs of the
   tree-level minimum as a new starting point for PyMinuit, & then compares the
   point to which PyMinuit rolled from this new starting point to the input
   minimum (as basins of attraction can sometimes move so far with loop
   corrections that a tree-level CCB minimum can get caught in the input
   minimum's basin of attraction).
 - Default Python now also now restricts PyMinuit to a hypercube of each field
   being within a hundred times the scale of the SLHAfile in magnitude (rather
   than a thousand).

 * 8th July 2013: version 1.0.4
 - Updated citation text, as Vevacious has now got a manual on the arXiv!
   [arXiv:1307.1477 (hep-ph)]
 - Updated README

 * 28th June 2013: version 1.0.3
 - Fixed default Vevacious.py to correctly handle PyMinuit exceptions, and to
   correctly warn when the input VEVs seem to be giving a saddle point.

 * 26th June 2013: version 1.0.2
 - Fixed SARAH-generated example model files to no longer have flavor issue, so
   now they are as they would be if generated by the current version of SARAH4.

 * 26th June 2013: version 1.0.1
 - Fixed "pure SLHA" example model files.
 - Updated NMSSM_wrong_neutral example point.
 
 * 17th June 2013: version 1.0.0
 - Release version!
 - Added consistency check that all used SLHA BLOCK scales must be the same,
   otherwise the program prints an error message and returns EXIT_FAILURE.
 - Updated examples.
   
 * 5th June 2013: version 0.3.1
 - Added VevaciousRunner::findTreeLevelExtrema(...) as a synonym for
   VevaciousRunner::overwriteTreeLevelExtrema(...).
 - Changed verdict of results from "unstable" for tunneling times below the
   threshold to "short-lived" and from "metastable" to "long-lived" for those
   over the threshold.
 - Added example model files in Vevacious/MSSM/ that use purely the SLHA1
   conventions without requiring the extra HMIX data lines that are assumed by
   SARAH.

 * 17th May 2013: version 0.3.0
 - Vevacious now prints its version and the documentation citation in the
   results file and also SLHA block.
 - Added NMSSM example model files and a single, metastable example spectrum.

 * 14th May 2013: version 0.2.9
 - First version publicly available, but not officially released yet.
 - Slight changes to code calling CosmoTransitions to no longer use
   quickTunneling = True, on advice from the author of CosmoTransitions, Max
   Wainwright.
 - saddle_nudges input argument added, to give a list of floating point numbers
   to be used as the set of nudge sizes for the displacement of MINUIT off
   saddle points.
 - Default Vevacious.py cleaned up a bit to make it easier to substitute in the
   tree-level potential as the function to be used for the analysis. PyMinuit
   also has limits of a thousand times the energy scale for the VEVs now, to
   avoid it rolling off to infinity.
   
 * 26th April 2013: version 0.2.8
 - Publicly available on GitHub, but still not officially released!
 - added optional bool argument to VevaciousRunner::appendResultsToSlha
   which inserts '#' after the doubles of the VEVACIOUSRESULTS block so that
   SSP can read it without problems, & added argument option to Vevacious.exe
   so that this can be used by default.
   
 * 17th April 2013: version 0.2.7
 - Still not even released!
 - Oops, fixed a bug from an undocumented change which was designed to prevent
   CosmoTransitions from trying to tunnel through an energy barrier that is
   smaller than its resolution, which was incorrectly deciding based on whether
   or not the 1st step was absolutely positive, rather than relative to the
   input VEVs.
 - Tried to add use of BOL::WaitingOnSubprocessExecutor to run Python with
   intervention to kill the process if it took too long, so that it could be
   run again without allowing path deformation, but it doesn't work for arcane
   path-inheritance-for-subprocess issues.
   
 * 17th April 2013: version 0.2.6
 - Still not even released!
 - Moved appending SLHA blocks to
   VevaciousRunner::appendResultsToSlha( std::string const& ) const function
   out of Vevacious.cpp, which now calls this function.
   
 * 16th April 2013: version 0.2.5
 - Still not even released!
 - Will now format Mathematica-style 123.4e5 into 123.4E+5 so that HOM4PS2 is
   happy, & other such numbers where 'e' or 'E' is not followed by '+' or '-'.
   
 * 2nd April 2013: version 0.2.4
 - Still not even released!
 - Fixed bug where exceptions from PyMinuit were not being caught correctly,
   leading to no results being printed from a parameter point that has a
   PyMinuit starting point that causes an exception.
 - Changed behavior such that result file is removed before doing any
   calculation, so that an error does not leave an old result file.

 * 29th March 2013: version 0.2.3
 - Still not even released!
 - Fixed bug where PotentialMinimizer::prepareLoopCorrections (and thus
   VevaciousRunner::prepareLoopCorrections) would only work for the 1st SLHA
   file given, since I had forgotten to clear the list of mass correction
   function names.
 - Tidied up by renaming remaining references to tree-level potentials into
   references to polynomial parts of the potential.

 * 27th March 2013: version 0.2.2
 - Still not even released!
 - Updated READMEs & example files.
 
 * 25th March 2013: version 0.2.1
 - Still not even released!
 - Before any .py file is written, Vevacious deletes any pre-existing .pyc
   version of the file, allowing multiple parameter points to be run per
   second.
 - In the default Vevacious.py, the minimum found by MINUIT after being rolled
   from the input VEVs, which is used as the actual input minimum, is made
   deeper by twice the MINUIT error, so that numerically degenerate minima are
   not considered as global minimum candidates.

 * 22nd March 2013: version 0.2.0
 - Still not even released!
 - Now tries to find blocks with prefixes for "tree-level" and "1-loop-level"
   for corresponding parts to be consistent with the renormalization scheme
   used by SPheno. By default the "tree-level" prefix is "TREE" and the
   "1-loop-level" prefix is "LOOP", but these can be changed in the model file
   by using the block_prefixes XML element to set the 'tree' and 'loop'
   attributes (if the element is included in the file, attributes not
   explicitly written are assumed to be empty strings). The tadpole equations
   & the mass-squared matrices use "tree-level values" if available, and the
   polynomial part of the potential uses "1-loop-level" values if available.
   E.g. for the mu parameter, given in HMIX[1]: in the tadpole equations,
   TREEHMIX[1] is used if written in the SLHA file, & if it is not found,
   HMIX[1] is used; in the polynomial part of the potential, LOOPHMIX[1] is
   used if found in the SLHA file, HMIX[1] otherwise.

 * 15th March 2013: version 0.1.5
 - Still not even released!
 - Now the parameter-dependent Python will be properly written even if HOM4PS2
   has not been run (was not correctly writing FunctionFromDictionary and other
   functions that require knowledge of the internal VEV names).
   
 * 14th March 2013: version 0.1.4
 - Still not even released!
 - Fixed bug where ./Vevacious.py wasn't being written, instead
   Vevacious.py.stau_VEVs_MSSM was always being written (but that was just a
   name, its content was completely independent of the model).
 - Now the name of the Python main program to run (by default, ./Vevacious.py)
   is given by the python_main element of the initialization XML file. If there
   is no file with that name in existence, the default ./Vevacious.py will be
   written and run.

 * 12th March 2013: version 0.1.3
 - Still not even released!
 - Now allows for a tolerance on the imaginary parts of the HOM4PS2 solutions
   to be considered valid real solutions.
 - Now appends the warnings to the SLHA file as the VEVACIOUSWARNINGS block.

 * 11th March 2013: version 0.1.2
 - Still not even released!
 - Default Vevacious.py no longer tries extrema that are not as deep as the
   input but happen to have a depth difference within the MINUIT error.
 - Default Vevacious.py now records warnings in the result file.
 - Default Vevacious.py now copes with PyMinuit failing to find a minimum with
   sufficiently positive Hessian eigenvalues.

 * 10th March 2013: version 0.1.1
 - Still not even released! Credits updated.
 - Rewritten such that the tree-level extrema are written to a Python file,
   and the parameter-point dependent stuff (functions, VEV name maps) are
   written to another file, and the minimization is performed by running a
   parameter-point-independent Python program, either supplied by the user or
   a default program is written. (The intent is that the user can modify the
   main Python code without having to recompile any C++.)

 * 2nd January 2013: version 0.0.1
 - Not even released! Credits updated.
 - Fixed missing factor of (1/(16 pi^2)) from loop corrections.


The C++ files of Vevacious are:

 <> headers in Vevacious/include/:
 PotentialMinimizer.hpp
 SarahInterpreter.hpp
 SarahSlhaConverter.hpp
 TadpoleSolver.hpp
 VevaciousRunner.hpp
 VevRenamer.hpp

 <> source files in Vevacious/source/:
 PotentialMinimizer.cpp
 SarahInterpreter.cpp
 SarahSlhaConverter.cpp
 TadpoleSolver.cpp
 Vevacious.cpp
 VevaciousRunner.cpp
 VevRenamer.cpp

 <> and also:
 Vevacious/Makefile
 and README.Vevacious.txt which describes the package (copied as README.txt).
 The makefile creates a single executable in Vevacious/bin/:
 Vevacious.exe, the executable which performs the full set of steps to try to
 find the global minimum of the loop-corrected potential.
 
 <> Some example files are also given:
 * Vevacious/bin/VevaciousInitialization.xml which is a default initialization
 file, which uses paths that need to be changed.
 * Various files in Vevacious/MSSM/:
 ** MSSM_XYZ.vin is a model file for the MSSM using the VEVs specified by XYZ.
 For example, MSSM_RealHiggsAndStauVevs.vin is a model file for the MSSM
 allowing VEVs for the staus as well as the Higgs, with the assumption that
 they are all real.
 *** MSSM_JustNormalHiggsVevs.vin: just the normal vd and vu are allowed to be
 non-zero.
 *** MSSM_RealHiggsAndStauVevs.vin: in addition to the normal vd and vu, stau_L
 and stau_R are allowed to have non-zero VEVs.
 *** MSSM_RealHiggsAndStauAndStopVevs.vin: as MSSM_RealHiggsAndStauVevs.vin,
 but one color of stop_L and the same color of stop_R are allowed non-zero
 VEVs.
 ** pure_SLHA_XYZ.vin: versions of the above, but allowing pure SLHA2 input
 without requiring the additional information given by SPhenoMSSM.
 ** XYZ.slha.in is an input file for SPhenoMSSM (the SPheno made by SARAH for
 the model MSSM), XYZ.slha.out is an output spectrum from that input.
 *** XYZ = SPS1ap: the CMSSM parameter point SPS1a'
 *** XYZ = CMSSM_CCB: an example charge- and color-breaking global minimum
 CMSSM parameter point (corresponding to the best-fit point of arXiv:1204.4199)
 *** XYZ = NUHM1_CCB: an example charge- and color-breaking global minimum
 NUHM1 (Non-Universal Higgs Mass CMSSM with 1 additional mass for the Higgses,
 rather than 2) parameter point (corresponding to the "low" best-fit point of
 arXiv:1207.7315)
 * Various files in Vevacious/NMSSM/:
 As above, but for the NMSSM, though with only one example CNMSSM parameter
 point, corresponding to P1 of arXiv:0801.4321, which has a global minimum with
 the wrong Z mass, but still zero stau and stop VEVs.
 