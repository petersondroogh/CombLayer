/********************************************************************* 
  CombLayer : MCNP(X) Input builder
 
 * File:   singleItemBuild/singleItemVariables.cxx
 *
 * Copyright (c) 2004-2019 by Stuart Ansell
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>. 
 *
 ****************************************************************************/
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>
#include <list>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <algorithm>
#include <iterator>
#include <memory>

#include "Exception.h"
#include "FileReport.h"
#include "NameStack.h"
#include "RegMethod.h"
#include "GTKreport.h"
#include "OutputLog.h"
#include "support.h"
#include "MatrixBase.h"
#include "Matrix.h"
#include "Vec3D.h"
#include "Code.h"
#include "varList.h"
#include "FuncDataBase.h"
#include "variableSetup.h"
#include "singleItemVariables.h"

#include "CryoGenerator.h"
#include "BladeGenerator.h"
#include "RectPipeGenerator.h"
#include "TwinBaseGenerator.h"
#include "TwinGenerator.h"
#include "TwinFlatGenerator.h"
#include "ChopperGenerator.h"
#include "QuadrupoleGenerator.h"
#include "EPSeparatorGenerator.h"
#include "PreDipoleGenerator.h"
#include "DipoleChamberGenerator.h"

namespace setVariable
{


void
SingleItemVariables(FuncDataBase& Control)
  /*!
    Function to set the control variables and constants
    -- This version is for ESS ()
    \param Control :: Function data base to add constants too
  */
{
// -----------
// GLOBAL stuff
// -----------
  Control.addVariable("zero",0.0);     // Zero
  Control.addVariable("one",1.0);      // one

  // photon test
  Control.addVariable("TargetZStep",0.0);
  Control.addVariable("TargetRadius",5.0);
  Control.addVariable("TargetMat","Stainless304");

  Control.addVariable("ShieldInnerRadius",50.0);
  Control.addVariable("ShieldOuterRadius",100.0);
  Control.addVariable("ShieldMat","Stainless304");

  Control.addVariable("TubeYStep",10.0);
  Control.addVariable("TubeRadius",5.0);
  Control.addVariable("TubeLength",30.0);
  Control.addVariable("TubeDepth",10.0);
  Control.addVariable("TubeWidth",20.0);
  Control.addVariable("TubeHeight",20.0);
  Control.addVariable("TubeMat","Stainless304");
  
  Control.addVariable("TubeARadius",5.0);
  Control.addVariable("TubeALength",25.0);
  Control.addVariable("TubeBRadius",5.0);
  Control.addVariable("TubeBLength",5.0);
  Control.addVariable("TubeAMat","Stainless304");
  Control.addVariable("TubeBMat","Lead");
  
  Control.addVariable("CryoBOuterRadius",20.0);
  
  Control.addVariable("CryoBNLayers",0);
  Control.addVariable("CryoBLRadius0",1.0);
  Control.addVariable("CryoBLRadius1",3.0);
  Control.addVariable("CryoBLRadius2",5.0);
  Control.addVariable("CryoBLRadius3",8.0);
  Control.addVariable("CryoBLThick",0.5);
  Control.addVariable("CryoBLTemp",300.0);
  
  Control.addVariable("CryoBTopOffset",1.0);
  Control.addVariable("CryoBBaseOffset",1.0);
  Control.addVariable("CryoBCutTopAngle",10.0);
  Control.addVariable("CryoBCutBaseAngle",10.0);

  Control.addVariable("CryoBTopThick",10.0);
  Control.addVariable("CryoBBaseThick",10.0);

  Control.addVariable("CryoBApertureWidth",2.0);
  Control.addVariable("CryoBApertureHeight",2.0);

  Control.addVariable("CryoBMat","Aluminium");


  
  setVariable::CryoGenerator CryGen;
  CryGen.generateFridge(Control,"singleCryo",3.0,-10,4.5);

  setVariable::TwinGenerator TGen;
  TGen.generateChopper(Control,"singleTwinB",0.0,16.0,10.0);  

  setVariable::TwinFlatGenerator TCGen;
  TCGen.generateChopper(Control,"singleTwinC",0.0,16.0,10.0);  

  TCGen.generateChopper(Control,"singleTwinD",80.0,16.0,10.0);  

  setVariable::RectPipeGenerator PipeGen;
  PipeGen.generatePipe(Control,"singleBoxPipeA",0.0,80.0);
  PipeGen.generatePipe(Control,"singleBoxPipeB",0.0,80.0);

  setVariable::ChopperGenerator CGen;
  CGen.setMotorRadius(10.0);
  CGen.generateChopper(Control,"singleChopper",10.0,12.0,5.55);
  Control.addVariable("singleChopperMotorBodyLength",15.0);
  

  setVariable::BladeGenerator BGen;
  // Single Blade chopper
  BGen.setThick({0.2});
  BGen.addPhase({95,275},{30.0,30.0});
  BGen.generateBlades(Control,"singleBBladeTop",-2.0,22.5,35.0);

  // Single Blade chopper
  BGen.setThick({0.2});
  BGen.addPhase({95,275},{30.0,30.0});
  BGen.generateBlades(Control,"singleBBladeLow",2.0,22.5,35.0);

  setVariable::QuadrupoleGenerator QGen;
  QGen.generateQuad(Control,"Quad",0.0,10.0);

  setVariable::EPSeparatorGenerator EPGen;
  EPGen.generatePipe(Control,"EPSep",0.0);

  setVariable::PreDipoleGenerator PDGen;
  PDGen.generatePipe(Control,"PreDipole",0.0);

  setVariable::DipoleChamberGenerator DCGen;
  DCGen.generatePipe(Control,"DipoleChamber",0.0);
  
  
  return;
}

}  // NAMESPACE setVariable
