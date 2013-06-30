/*****************************************************************************\
 * README..txt                                                               *
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
 *      A full list of the files of HoTTMiLC is at the end of this file.     *
\*****************************************************************************/

 Now that the legalese preamble is out of the way, the description of the code
 and how to use it can begin!

 Vevacious is a program written in C++ to try to find the global minimum of a
 loop-corrected potential. It is designed to be used in conjuction with SARAH
 by Florian Staub and utilizes HOM4PS2 by Tsung-Lin Lee, Tien-Yien Li, and
 Chih-Hsiung Tsai, and PyMinuit by Jim Pivarski, implementing the MINUIT
 algorithms by Fred James.
 
 Vevacious is an implementation of the process of attempting to find the global
 minimum of a potential as decided upon by the collaboration of
 José Eliel Camargo (elielx@gmail.com),
 Ben O'Leary (benjamin.oleary@gmail.com),
 Werner Porod (porod@physik.uni-wuerzburg.de), and
 Florian Staub (florian.staub@googlemail.com). This README file was written by
 Ben O'Leary.
 
 Unfortunately, automating the full installation of the above codes is beyond
 my capabilities, so one has to follow the installation instructions given
 after the code description below (or wing it oneself).
 
 Vevacious runs in 2 stages: first it prepares the tree-level tadpole equation
 system and then solves it with HOM4PS2, and second, it uses the tree-level
 solutions as starting points for minimizing the 1-loop-level potential with
 PyMinuit.
 
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
 in a certain renormalization scheme)).
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
    (link last checked 2013-05-14). The installation is as simple as unpacking
    the gzipped tarball.
 2) Download HOM4PS2. The files are available at
    http://www.math.nsysu.edu.tw/~leetsung/works/HOM4PS_soft_files/HOM4PS_Linux.htm
    (link last checked 2013-05-14). The installation is as simple as unpacking
    the gzipped tarball.
 3) Ensure that Python is installed. PyMinuit requires at least version 2.4 or
    later. I shouldn't have to get into any specifics of how to install Python
    here... Internet search engines are your friends.
 4) Download and install PyMinuit. The instructions are available at
    http://code.google.com/p/pyminuit/wiki/HowToInstall
    (link last checked 2013-05-14). The installation involves downloading and
    compiling one of the C++ implementations of MINUIT (such as the standalone
    version 1.7.9 or the CERN ROOT version, both available from
    http://lcgapp.cern.ch/project/cls/work-packages/mathlibs/minuit/
    (link last checked 2013-05-14)), then running the PyMinuit setup script
    (giving the path where the C++ MINUIT code was _built_, *not* installed -
    it needs the .o object files rather than the .a library file...). The
    LD_LIBRARY_PATH and PYTHONPATH environment variables then need to be set.
 5) Download and compile the LesHouchesParserClasses (LHPC) C++ library. The
    files are available at
    http://www.hepforge.org/downloads/lhpc
    (link last checked 2013-05-14) or
    https://github.com/benoleary/LesHouchesParserClasses_CPP
    (link last checked 2013-05-14). The installation should just be
    downloading, unzipping, and then running make.
 6) Compile Vevacious. The Makefile should be edited to have the correct paths
    to the header files and library file of LHPC (these are
    /path/to/unzipped/LHPC/include/ and /path/to/unzipped/LHPC/lib/ by
    default). LHPC version 0.8.4 or higher is required.

 Once Vevacious is installed, it can be run as is. Various parameters for each
 execution can be specified: for the list, see the
 Vevacious/bin/VevaciousInitialization.xml file, where the various tags <XYZ>
 can be over-ridden by being given as arguments (e.g. --XYZ=ABC will over-ride
 whatever value is given between <XYZ> and </XYZ> in the .ini file). A
 different initialization file can be given with the option
 --input=other_filename as well (and direct arguments will over-ride any given
 in that file too).


CHANGELOG:
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
 Vevacious/bin/VevaciousInitialization.xml which is a default initialization
 file, which uses paths that need to be changed.
 Various files in Vevacious/MSSM/:
 Vevacious.in.XYZ is a model file for XYZ pregenerated from a beta version of
 SARAH4. In principle most VEVs should be taken as complex, which would require
 them to be written in the form vReal + i * vImag, but these examples, because
 of computational limitations, assume only real VEVs. The differences are which
 fields are allowed to have non-zero (real) VEVs:
 * XYZ = MSSM:
     just the neutral components of the Higgs doublets
 * XYZ = realStauVevs_MSSM:
     the neutral components of the Higgs doublets,
     stau_L,
     stau_R.
 * XYZ = realStauAndStopVevs_MSSM:
     the neutral components of the Higgs doublets,
     stau_L,
     stau_R,
     stop_L (1 color),
     stop_R (the same color as stop_L).

 
 
