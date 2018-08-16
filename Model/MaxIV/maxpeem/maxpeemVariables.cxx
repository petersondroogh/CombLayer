/********************************************************************* 
  CombLayer : MCNP(X) Input builder
 
 * File:   maxpeem/maxpeemVariables.cxx
 *
 * Copyright (c) 2004-2018 by Stuart Ansell
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
 * along with this program.  If not, see <http://www.gnu.org/licenses/>. 
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

#include "CFFlanges.h"
#include "PipeGenerator.h"
#include "SplitPipeGenerator.h"
#include "BellowGenerator.h"
#include "LeadPipeGenerator.h"
#include "CrossGenerator.h"
#include "GateValveGenerator.h"
#include "JawValveGenerator.h"
#include "PipeTubeGenerator.h"
#include "PortTubeGenerator.h"
#include "PortItemGenerator.h"
#include "VacBoxGenerator.h"
#include "FlangeMountGenerator.h"
#include "MirrorGenerator.h"
#include "CollGenerator.h"
#include "SqrFMaskGenerator.h"
#include "PortChicaneGenerator.h"
#include "LeadBoxGenerator.h"

namespace setVariable
{

namespace maxpeemVar
{

void
heatDumpVariables(FuncDataBase& Control,const std::string& frontKey)
  /*!
    Build the heat dump variables
    \param Control :: Database
    \param frontKey :: prename
   */
{
  ELog::RegMethod RegA("balderVariables","heatDumpVariables");

  setVariable::PortTubeGenerator PTubeGen;
  setVariable::PortItemGenerator PItemGen;
  setVariable::FlangeMountGenerator FlangeGen;

  PTubeGen.setMat("Stainless304");
  PTubeGen.setCF<CF150>();
  PTubeGen.setPortLength(2.5,2.5);
  
  PTubeGen.generateCFTube<CF150>(Control,frontKey+"HeatBox",40.0,20.0);
  Control.addVariable(frontKey+"HeatBoxZAngle",90);
  Control.addVariable(frontKey+"HeatBoxNPorts",2);

  // beam ports
  PItemGen.setCF<setVariable::CF40>(5.0);
  PItemGen.setPlate(0.0,"Void");  

  FlangeGen.setCF<setVariable::CF40>();
  FlangeGen.setBlade(5.0,10.0,1.0,0.0,"Tungsten",0);     // W / H / T

  const Geometry::Vec3D ZVec(0,0,1);
  const std::string heatName=frontKey+"HeatBoxPort";
  const std::string hName=frontKey+"HeatDumpFlange";
  PItemGen.generatePort(Control,heatName+"0",Geometry::Vec3D(0,0,0),ZVec);
  PItemGen.generatePort(Control,heatName+"1",Geometry::Vec3D(0,0,0),-ZVec);


  const std::string hDump(frontKey+"HeatDump");
  Control.addVariable(hDump+"Height",10.0);
  Control.addVariable(hDump+"Width",3.0);
  Control.addVariable(hDump+"Thick",8.0);
  Control.addVariable(hDump+"CutHeight",10.0);
  Control.addVariable(hDump+"CutDepth",0.0);
  Control.addVariable(hDump+"Mat","Tungsten");
  return;
}
  
void
frontEndVariables(FuncDataBase& Control,
		  const std::string& frontKey)
/*!
    Set the variables for the mono
    \param Control :: DataBase to use
    \param frontKey :: name before part names
  */
{
  ELog::RegMethod RegA("balderVariables[F]","frontEndVariables");

  setVariable::BellowGenerator BellowGen;
  setVariable::PipeGenerator PipeGen;
  setVariable::CrossGenerator CrossGen;
  setVariable::PipeTubeGenerator SimpleTubeGen;
  setVariable::VacBoxGenerator VBoxGen;
  setVariable::SqrFMaskGenerator CollGen;
  setVariable::PortTubeGenerator PTubeGen;
  setVariable::PortItemGenerator PItemGen;
  setVariable::FlangeMountGenerator FlangeGen;
  
  PipeGen.setWindow(-2.0,0.0);   // no window
  PipeGen.setMat("Stainless304");
  
  VBoxGen.setMat("Stainless304");
  VBoxGen.setWallThick(1.0);
  VBoxGen.setCF<CF40>();
  VBoxGen.setPortLength(5.0,5.0); // La/Lb
  // ystep/width/height/depth/length
  VBoxGen.generateBox(Control,frontKey+"WigglerBox",
		      115.0,30.0,15.0,15.0,210.0);

  // Wiggler
  Control.addVariable(frontKey+"WigglerLength",200.0);
  Control.addVariable(frontKey+"WigglerBlockWidth",8.0);
  Control.addVariable(frontKey+"WigglerBlockHeight",8.0);
  Control.addVariable(frontKey+"WigglerBlockDepth",8.0);
  Control.addVariable(frontKey+"WigglerBlockHGap",0.2);
  Control.addVariable(frontKey+"WigglerBlockVGap",0.96);

  Control.addVariable(frontKey+"WigglerBlockVCorner",1.0);
  Control.addVariable(frontKey+"WigglerBlockHCorner",2.0);

  
  Control.addVariable(frontKey+"WigglerVoidMat",0);
  Control.addVariable(frontKey+"WigglerBlockMat","Iron_10H2O");

  Control.addVariable(frontKey+"ECutDiskYStep",2.0);
  Control.addVariable(frontKey+"ECutDiskLength",0.1);
  Control.addVariable(frontKey+"ECutDiskRadius",0.11);
  Control.addVariable(frontKey+"ECutDiskDefMat","H2Gas#0.1");

  PipeGen.setCF<CF40>();
  PipeGen.generatePipe(Control,frontKey+"DipolePipe",0,444.50);

  BellowGen.setCF<setVariable::CF40>();
  BellowGen.setBFlangeCF<setVariable::CF63>();
  BellowGen.generateBellow(Control,frontKey+"BellowA",0,10.0);

  // collimator block
  CollGen.setCF<CF63>();
  CollGen.setBFlangeCF<CF40>();
  CollGen.setFrontGap(3.99,1.97);  //1033.8
  CollGen.setBackGap(0.71,0.71);
  CollGen.setMinSize(10.2,0.71,0.71);
  CollGen.generateColl(Control,frontKey+"CollA",0.0,15.0);


  BellowGen.setCF<setVariable::CF40>();
  BellowGen.generateBellow(Control,frontKey+"BellowB",0,10.0);

  // flange if possible
  CrossGen.setPlates(0.5,2.0,2.0);  // wall/Top/base
  CrossGen.setTotalPorts(10.0,10.0);     // len of ports (after main)
  CrossGen.setMat("Stainless304");

  // height/depth
  CrossGen.generateDoubleCF<setVariable::CF40,setVariable::CF100>
    (Control,frontKey+"IonPA",0.0,26.6,26.6);

  BellowGen.setCF<setVariable::CF40>();
  BellowGen.generateBellow(Control,frontKey+"BellowC",0,10.0);

  PipeGen.setCF<CF40>();
  PipeGen.generatePipe(Control,frontKey+"HeatPipe",0,115.0);

  PipeGen.setCF<CF40>();
  PipeGen.generatePipe(Control,frontKey+"HeatPipe",0,115.0);

  heatDumpVariables(Control,frontKey);
  return;
}

}  // NAMESPACE maxpeemVar

  
void
MAXPEEMvariables(FuncDataBase& Control)
  /*!
    Function to set the control variables and constants
    -- This version is for Photon Moderator
    \param Control :: Function data base to add constants too
  */
{
  ELog::RegMethod RegA("maxpeemVariables[F]","maxpeemVariables");

  Control.addVariable("sdefType","Wiggler");

  maxpeemVar::frontEndVariables(Control,"MaxPeemFrontBeam");  

  return;
}

}  // NAMESPACE setVariable
