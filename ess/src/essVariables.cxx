/********************************************************************* 
  CombLayer : MNCPX Input builder
 
 * File:   essBuild/essVariables.cxx
*
 * Copyright (c) 2004-2013 by Stuart Ansell
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
#include <iterator>
#include <boost/array.hpp>
#include <boost/shared_ptr.hpp>

#include "Exception.h"
#include "FileReport.h"
#include "NameStack.h"
#include "RegMethod.h"
#include "GTKreport.h"
#include "OutputLog.h"
#include "support.h"
#include "stringCombine.h"
#include "MatrixBase.h"
#include "Matrix.h"
#include "Vec3D.h"
#include "Code.h"
#include "varList.h"
#include "FuncDataBase.h"
#include "variableSetup.h"

namespace setVariable
{
  void EssBeamLinesVariables(FuncDataBase&);
  void EssConicModerator(FuncDataBase&);
  void ESSWheel(FuncDataBase&);
  void ESSLayerMod(FuncDataBase&);

void ESSLayerMod(FuncDataBase& Control)
  /*
    Set the variables for the lower moderators
    \param Control :: DataBase to put variables
   */
{
  //  Control.addVariable("LowModXStep", 10.0);  
  //  Control.addVariable("LowModYStep", 10.0);  
  Control.addVariable("LowModZStep",-15.85);
  Control.addVariable("LowModXYangle",125.15); 
  Control.addVariable("LowModZangle",0.0);
  Control.addVariable("LowModRadius",8.0);
  Control.addVariable("LowModHeight", 8.0);
  Control.addVariable("LowModMat", "M01001");
  Control.addVariable("LowModTemp",20.0);
  Control.addVariable("LowModNLayers",7);
  // al layer
  Control.addVariable("LowModHGap1",0.3);
  Control.addVariable("LowModRadGap1",0.3);
  Control.addVariable("LowModMaterial1", 13061);  // Al materk
  Control.addVariable("LowModTemp1",20.0);  
  // Vac gap
  Control.addVariable("LowModHGap2",0.5);
  Control.addVariable("LowModRadGap2",0.5);
  Control.addVariable("LowModMaterial2",0); 
  // Next Al layer
  Control.addVariable("LowModHGap3",0.2);
  Control.addVariable("LowModRadGap3",0.2);
  Control.addVariable("LowModMaterial3", "Aluminium"); 
  //Control.addVariable("LowModTemp3",400);//77.0);   kbat
  // He Layer
  Control.addVariable("LowModHGap4",0.2);
  Control.addVariable("LowModRadGap4",0.2);
  Control.addVariable("LowModMaterial4",0); 
  // Outer Layer
  Control.addVariable("LowModHGap5",0.2);
  Control.addVariable("LowModRadGap5",0.2);
  Control.addVariable("LowModMaterial5","Aluminium"); 
  //  Control.addVariable("LowModTemp5",300.0); 
  // Clearance
  Control.addVariable("LowModHGap6",0.2);
  Control.addVariable("LowModRadGap6",0.2);
  Control.addVariable("LowModMaterial6",0); 
  return;
}



void
EssWheel(FuncDataBase& Control)
  /*!
    Variables that are used for the segmented wheel
    \param Control :: Segment variables
   */
{
  // WHEEL SHAFT

  Control.addVariable("SegWheelShaftNLayers",6);
  Control.addVariable("SegWheelShaftTopHeight",400.0);
  Control.addVariable("SegWheelShaftBaseHeight",-13.5);
  // Control.addVariable("SegWheelShaftHeight",435.0);

  // Control.addVariable("SegWheelShaftRadius",25.0);
  // see WheelInnerRadius below (28cm): it's the external joint radius
  Control.addVariable("SegWheelShaftJointThick",3.0);
  Control.addVariable("SegWheelShaftCoolThick",4.0);
  Control.addVariable("SegWheelShaftCladThick",3.0);
  Control.addVariable("SegWheelShaftVoidThick",1.0);
  Control.addVariable("SegWheelShaftSupportRadius",40.0);
  Control.addVariable("SegWheelShaftSupportThick",2.5);
  Control.addVariable("SegWheelShaftBaseThick",3.0);
  Control.addVariable("SegWheelShaftBaseFootThick",13.5);

  Control.addVariable("SegWheelCladShaftMat","Stainless304");
  Control.addVariable("SegWheelCoolingShaftMatInt",26316);
  Control.addVariable("SegWheelCoolingShaftMatExt",0);
 
  // TARGET

  Control.addVariable("SegWheelXStep",0.0);  
  Control.addVariable("SegWheelYStep",115.0);  
  Control.addVariable("SegWheelZStep",0.0);
  Control.addVariable("SegWheelXYangle",0.0); 
  Control.addVariable("SegWheelZangle",0.0);
  //
  Control.addVariable("SegWheelTargetHeight",8.0);
  Control.addVariable("SegWheelTargetSectorOffsetX",0.0);  
  Control.addVariable("SegWheelTargetSectorOffsetY",142.6476749356);
  Control.addVariable("SegWheelTargetSectorOffsetZ",0.0);
  Control.addVariable("SegWheelTargetSectorAngleXY",1.1234180949);
  Control.addVariable("SegWheelTargetSectorAngleZ",0.0);
  Control.addVariable("SegWheelTargetSectorApertureXY",2.2468361899);
  Control.addVariable("SegWheelTargetSectorApertureZ",0.0);
  Control.addVariable("SegWheelTargetSectorNumber",33);
  //
  Control.addVariable("SegWheelCoolantThickOut",0.15);
  Control.addVariable("SegWheelCoolantThickIn",0.65);
  //
  Control.addVariable("SegWheelCaseThickZ",0.5);
  Control.addVariable("SegWheelCaseThickX",0.1);
  //
  Control.addVariable("SegWheelVoidThick",2.15);  //distance target-ion tube height

  Control.addVariable("SegWheelInnerRadius",28.0);
  Control.addVariable("SegWheelCoolantRadiusOut",124.8);
  Control.addVariable("SegWheelCoolantRadiusIn",114.5);
  Control.addVariable("SegWheelCaseRadius",125.0);
  Control.addVariable("SegWheelVoidRadius",126.0);

  // Material types 1:2:3
  Control.addVariable("SegWheelWMat","Tungsten");
  Control.addVariable("SegWheelSteelMat","Stainless304");
  Control.addVariable("SegWheelHeMat","helium");

  Control.addVariable("SegWheelInnerMat","Stainless304");

  Control.addVariable("SegWheelNLayers",25);

  Control.addVariable("SegWheelRadius1",84.27);
  Control.addVariable("SegWheelMatTYPE1",1);

  Control.addVariable("SegWheelRadius2",85.65);
  Control.addVariable("SegWheelMatTYPE2",2);

  Control.addVariable("SegWheelRadius3",97.02);
  Control.addVariable("SegWheelMatTYPE3",26317);

  Control.addVariable("SegWheelRadius4",97.52);
  Control.addVariable("SegWheelMatTYPE4",2);

  Control.addVariable("SegWheelRadius5",102.37); 
  Control.addVariable("SegWheelMatTYPE5",26317);

  Control.addVariable("SegWheelRadius6",102.77);
  Control.addVariable("SegWheelMatTYPE6",2);

  Control.addVariable("SegWheelRadius7",106.97);
  Control.addVariable("SegWheelMatTYPE7",26317);

  Control.addVariable("SegWheelRadius8",107.29);
  Control.addVariable("SegWheelMatTYPE8",2);

  Control.addVariable("SegWheelRadius9",110.49);
  Control.addVariable("SegWheelMatTYPE9",26317);

  Control.addVariable("SegWheelRadius10",110.78);
  Control.addVariable("SegWheelMatTYPE10",2); 

  Control.addVariable("SegWheelRadius11",112.78);
  Control.addVariable("SegWheelMatTYPE11",26317); 

  Control.addVariable("SegWheelRadius12",113.05);
  Control.addVariable("SegWheelMatTYPE12",2); 

  Control.addVariable("SegWheelRadius13",114.65);
  Control.addVariable("SegWheelMatTYPE13",26317); 

  Control.addVariable("SegWheelRadius14",114.9);
  Control.addVariable("SegWheelMatTYPE14",2); 

  // Control.addVariable("SegWheelRadius15",114.9);
  // Control.addVariable("SegWheelMatTYPE15",1); 

  Control.addVariable("SegWheelRadius15",116.5);
  Control.addVariable("SegWheelMatTYPE15",26317); 

  Control.addVariable("SegWheelRadius16",116.75);
  Control.addVariable("SegWheelMatTYPE16",2); 

  Control.addVariable("SegWheelRadius17",118.35);
  Control.addVariable("SegWheelMatTYPE17",26317); 

  Control.addVariable("SegWheelRadius18",118.6);
  Control.addVariable("SegWheelMatTYPE18",2); 

  Control.addVariable("SegWheelRadius19",120.0);
  Control.addVariable("SegWheelMatTYPE19",26317); 

  Control.addVariable("SegWheelRadius20",120.25);
  Control.addVariable("SegWheelMatTYPE20",2); 

  Control.addVariable("SegWheelRadius21",121.65);
  Control.addVariable("SegWheelMatTYPE21",26317); 

  Control.addVariable("SegWheelRadius22",121.9);
  Control.addVariable("SegWheelMatTYPE22",2); 

  Control.addVariable("SegWheelRadius23",123.1);
  Control.addVariable("SegWheelMatTYPE23",26317); 

  Control.addVariable("SegWheelRadius24",123.35);
  Control.addVariable("SegWheelMatTYPE24",2); 

  Control.addVariable("SegWheelRadius25",124.55);
  Control.addVariable("SegWheelMatTYPE25",26317); 
  return;
}


void
EssVariables(FuncDataBase& Control)
  /*!
    Function to set the control variables and constants
    -- This version is for ESS ()
    \param Control :: Function data base to add constants too
  */
{
// -----------
// GLOBAL stuff
// -----------

  const double TopModHeight = 13 ; // para H content
  const double ModZStep = 15.85+2.12+0.03; // target-moderator distance
  const double BeRefXStep = 5.3729961E+00; // as in TSM130318
  const double BeRefYStep = 8.4339145E+00; // as in TSM130318
  const double TopFlightWidth = 26;
  const double TopPreABBlockHeight = 17.6;
  const double LowPreABBlockHeight = TopPreABBlockHeight; //16.2;
  const double PreABBlockLength = 11;

  // Golden rule - there must be these 2 variables. Changin them allows to scale the model.
  Control.addVariable("zero",0.0);     // Zero
  Control.addVariable("one",1.0);      // one

  // water cylinder
  Control.addVariable("LowWaterDiscXStep", 0.0);
  Control.addVariable("LowWaterDiscYStep", 0.0);  
  Control.addVariable("LowWaterDiscZStep", 0.2);
  Control.addVariable("LowWaterDiscXYangle",0.0); 
  Control.addVariable("LowWaterDiscZangle", 0.0);
  Control.addVariable("LowWaterDiscRadius",33.0);
  Control.addVariable("LowWaterDiscHeight", 3.3); // Water height=3 + wall height = 0.3
  Control.addVariable("LowWaterDiscWallThick", 0.3); // in addition to Be
  Control.addVariable("LowWaterDiscRefMat",    1011); 
  Control.addVariable("LowWaterDiscRefMat1",    0); 
  Control.addVariable("LowWaterDiscWidth",    -1); 
  Control.addVariable("LowWaterDiscLength",    -1); 
  Control.addVariable("LowWaterDiscWallMat",   13000);
  Control.addVariable("LowWaterDiscVoidThick", 0.0);

  Control.addVariable("TopWaterDiscXStep", 0.0);
  Control.addVariable("TopWaterDiscYStep", 0.0);  
  Control.addVariable("TopWaterDiscZStep", 0.2);
  Control.addVariable("TopWaterDiscXYangle",0.0); 
  Control.addVariable("TopWaterDiscZangle", 0.0);
  Control.addVariable("TopWaterDiscRadius",33.0);
  Control.addVariable("TopWaterDiscHeight", 3.3); // Water height=3 + wall height = 0.3
  Control.addVariable("TopWaterDiscWallThick", 0.3); // in addition to Be
  Control.addVariable("TopWaterDiscRefMat",    1011); 
  Control.addVariable("TopWaterDiscRefMat1",    0); 
  Control.addVariable("TopWaterDiscWidth",    -1); 
  Control.addVariable("TopWaterDiscLength",    -1); 
  Control.addVariable("TopWaterDiscWallMat",   13000);
  Control.addVariable("TopWaterDiscVoidThick", 0.0);
  Control.addVariable("TopWaterDiscVoidCellMat", 0);
  Control.addVariable("TopWaterDiscVoidCellHeight", 0.0);

  // onion cooling
  Control.addVariable("OnionXStep", 0);
  Control.addVariable("OnionYStep", 0);
  Control.addVariable("OnionZStep", 0.0);
  Control.addVariable("OnionXYangle",0.0); 
  Control.addVariable("OnionZangle",0.0);
  Control.addVariable("OnionHeight", 3); // should be same as moderator height
  Control.addVariable("OnionWallThick", 0.3); // in addition to Be
  Control.addVariable("OnionWallMat",   13001);
  Control.addVariable("OnionNRings", 2);
  Control.addVariable("OnionRadius1", 4);
  Control.addVariable("OnionGateWidth1", 1);
  Control.addVariable("OnionGateLength1", 2);
  Control.addVariable("OnionRadius2", 8);
  Control.addVariable("OnionGateWidth2", 2);
  Control.addVariable("OnionGateLength2", 2);

  Control.addVariable("OnionPreXStep", 0);
  Control.addVariable("OnionPreYStep", 0);
  Control.addVariable("OnionPreZStep", 0.0);
  Control.addVariable("OnionPreXYangle",0.0); 
  Control.addVariable("OnionPreZangle",0.0);
  Control.addVariable("OnionPreHeight", 10); // should be same as moderator+top/bottom premoderator height
  Control.addVariable("OnionPreWallThick", 0.3); // in addition to Be
  Control.addVariable("OnionPreWallMat",   13001);
  Control.addVariable("OnionPreNRings", 2);
  Control.addVariable("OnionPreRadius1", 4);
  Control.addVariable("OnionPreGateWidth1", 1);
  Control.addVariable("OnionPreGateLength1", 2);
  Control.addVariable("OnionPreRadius2", 8);
  Control.addVariable("OnionPreGateWidth2", 2);
  Control.addVariable("OnionPreGateLength2", 2);
  
  
  // low mod supply pipe
  Control.addVariable("LSupplyNSegIn",2);
  bool flat = !false;
  // Central point:
  if (flat) {
  // as in Flat
    Control.addVariable("LSupplyPPt0",Geometry::Vec3D(0, 0.0, 0.0));
    Control.addVariable("LSupplyPPt1",Geometry::Vec3D(0, 19,  0.0));
    Control.addVariable("LSupplyPPt2",Geometry::Vec3D(0, 19, -64));
    // Central point [Top]: !!! redefine coordinates
    Control.addVariable("LSupplyTopPPt0",Geometry::Vec3D(0,-1.0,5.0));
    Control.addVariable("LSupplyTopPPt1",Geometry::Vec3D(0,-19.25,5.0));
    Control.addVariable("LSupplyTopPPt2",Geometry::Vec3D(3.005,-19.25,64.930));
  } else {
    Control.addVariable("LSupplyPPt0",Geometry::Vec3D(0,-1.0,0.0));
    Control.addVariable("LSupplyPPt1",Geometry::Vec3D(0,-19.25,0.0));
    Control.addVariable("LSupplyPPt2",Geometry::Vec3D(3.005,-19.25,64.930));
  // Central point [Top]:
    Control.addVariable("LSupplyTopPPt0",Geometry::Vec3D(0,-1.0,5.0));
    Control.addVariable("LSupplyTopPPt1",Geometry::Vec3D(0,-19.25,5.0));
    Control.addVariable("LSupplyTopPPt2",Geometry::Vec3D(3.005,-19.25,64.930));
  }

  Control.addVariable("LSupplyNRadii",9);
  Control.addVariable("LSupplyRadius0",1.5);
  Control.addVariable("LSupplyRadius1",1.7);
  Control.addVariable("LSupplyRadius2",1.8);
  Control.addVariable("LSupplyRadius3",2.3);
  Control.addVariable("LSupplyRadius4",2.5);
  Control.addVariable("LSupplyRadius5",2.7);
  Control.addVariable("LSupplyRadius6",2.9);
  Control.addVariable("LSupplyRadius7",3.5);
  Control.addVariable("LSupplyRadius8",3.7);

  Control.addVariable("LSupplyMat0","ParaH2");
  Control.addVariable("LSupplyMat1","Aluminium");
  Control.addVariable("LSupplyMat2","Aluminium");
  Control.addVariable("LSupplyMat3","Void");
  Control.addVariable("LSupplyMat4","Aluminium");
  Control.addVariable("LSupplyMat5","Void");
  Control.addVariable("LSupplyMat6","Aluminium");
  Control.addVariable("LSupplyMat7","Void");
  Control.addVariable("LSupplyMat8","Aluminium");

  Control.addVariable("LSupplyTemp0",25.0);
  Control.addVariable("LSupplyTemp1",25.0);
  Control.addVariable("LSupplyTemp2",25.0);
  Control.addVariable("LSupplyTemp3",0.0);

  /*  Control.addVariable("LSupplyActive0",3);
  Control.addVariable("LSupplyActive1",7);
  Control.addVariable("LSupplyActive2",31);
  Control.addVariable("LSupplyActive3",127);
  Control.addVariable("LSupplyActive4",511);
  Control.addVariable("LSupplyActive5",255);
  Control.addVariable("LSupplyActive6",255);*/




  // low mod return pipe

  Control.addVariable("LReturnNSegIn",1);
  Control.addVariable("LReturnPPt0",Geometry::Vec3D(0,0,0));
  Control.addVariable("LReturnPPt1",Geometry::Vec3D(0,30,0));


  Control.addVariable("LReturnInRadius",1.5);
  Control.addVariable("LReturnInAlRadius",1.7);
  Control.addVariable("LReturnMidAlRadius",1.8);
  Control.addVariable("LReturnVoidRadius",2.3);
  Control.addVariable("LReturnOutAlRadius",2.5);

  Control.addVariable("LReturnInMat",1001);
  Control.addVariable("LReturnInAlMat", 13060);
  Control.addVariable("LReturnMidAlMat", 13060);
  Control.addVariable("LReturnVoidMat",0);
  Control.addVariable("LReturnOutAlMat", 13060);


  // top mod supply pipe
  Control.addVariable("TSupplyNSegIn",2);
  if (flat) {
    Control.addVariable("TSupplyPPt0",Geometry::Vec3D(0, -1, 0.0));
    Control.addVariable("TSupplyPPt1",Geometry::Vec3D(0, 19, 0.0));
    Control.addVariable("TSupplyPPt2",Geometry::Vec3D(0, 19, -190));
  } else { // tdr mod, vertical extraction
    // Central point:
    Control.addVariable("TSupplyPPt0",Geometry::Vec3D(0,-1.0,0.0));
    Control.addVariable("TSupplyPPt1",Geometry::Vec3D(0,-19.25,0.0));
    Control.addVariable("TSupplyPPt2",Geometry::Vec3D(3.005,-19.25,64.930));
  // Central point [Top]:
    Control.addVariable("TSupplyTopPPt0",Geometry::Vec3D(0,-1.0,5.0));
    Control.addVariable("TSupplyTopPPt1",Geometry::Vec3D(0,-19.25,5.0));
    Control.addVariable("TSupplyTopPPt2",Geometry::Vec3D(3.005,-19.25,64.930));
  }
  Control.addVariable("TSupplyNRadii",9);
  Control.addVariable("TSupplyRadius0",1.5);
  Control.addVariable("TSupplyRadius1",1.7);
  Control.addVariable("TSupplyRadius2",1.8);
  Control.addVariable("TSupplyRadius3",2.3);
  Control.addVariable("TSupplyRadius4",2.5);
  Control.addVariable("TSupplyRadius5",2.7);
  Control.addVariable("TSupplyRadius6",2.9);
  Control.addVariable("TSupplyRadius7",3.5);
  Control.addVariable("TSupplyRadius8",3.7);

  Control.addVariable("TSupplyMat0","ParaH2");
  Control.addVariable("TSupplyMat1","Aluminium");
  Control.addVariable("TSupplyMat2","Aluminium");
  Control.addVariable("TSupplyMat3","Void");
  Control.addVariable("TSupplyMat4","Aluminium");
  Control.addVariable("TSupplyMat5","Void");
  Control.addVariable("TSupplyMat6","Aluminium");
  Control.addVariable("TSupplyMat7","Void");
  Control.addVariable("TSupplyMat8","Aluminium");

  Control.addVariable("TSupplyTemp0",25.0);
  Control.addVariable("TSupplyTemp1",25.0);
  Control.addVariable("TSupplyTemp2",25.0);
  Control.addVariable("TSupplyTemp3",0.0);

  Control.addVariable("TSupplyActive0",3); // inner H2 and a layer of Al
  Control.addVariable("TSupplyActive1",7);
  Control.addVariable("TSupplyActive2",31);
  // flat geometry does not work with these lines uncommented: !!!
  //  Control.addVariable("TSupplyActive3",127);
  //  Control.addVariable("TSupplyActive4",511);
  //  Control.addVariable("TSupplyActive5",255);
  //  Control.addVariable("TSupplyActive6",255);


  // Control.addVariable("TSupplyNSegIn",2);
  // Control.addVariable("TSupplyPPt0",Geometry::Vec3D(0, 0, 0.0));
  // Control.addVariable("TSupplyPPt1",Geometry::Vec3D(0, -19.25, 0)); // height of the bend
  // // mod xyangle (125.5)+pipe xyangle (57.5) e.g 3.0052=65*SIN(-(125.15+57.5)/180*PI())
  // Control.addVariable("TSupplyPPt2",Geometry::Vec3D(3.005, -19.25, 64.930));

  // Control.addVariable("TSupplyInRadius",1.5);
  // Control.addVariable("TSupplyInAlRadius",1.7);
  // Control.addVariable("TSupplyMidAlRadius",1.8);
  // Control.addVariable("TSupplyVoidRadius",2.3);
  // Control.addVariable("TSupplyOutAlRadius",2.5);

  // Control.addVariable("TSupplyInMat", 1001); // bottom cold moderator - 'material in the tube'
  // Control.addVariable("TSupplyInAlMat",  13061);
  // Control.addVariable("TSupplyMidAlMat", 13061);
  // Control.addVariable("TSupplyVoidMat",0);
  // Control.addVariable("TSupplyOutAlMat", 13060);



  // top mod return pipe
  Control.addVariable("TReturnNSegIn",1);
  Control.addVariable("TReturnPPt0",Geometry::Vec3D(0,0,0));
  Control.addVariable("TReturnPPt1",Geometry::Vec3D(0,30,0));

  Control.addVariable("TReturnTopNSegIn",1);
  Control.addVariable("TReturnTopPPt0",Geometry::Vec3D(0,0,-5));
  Control.addVariable("TReturnTopPPt1",Geometry::Vec3D(0,30,-5));

  Control.addVariable("TReturnNRadii",8);
  Control.addVariable("TReturnRadius0",1.7);
  Control.addVariable("TReturnRadius1",1.8);
  Control.addVariable("TReturnRadius2",2.3);
  Control.addVariable("TReturnRadius3",2.5);
  Control.addVariable("TReturnRadius4",2.7);
  Control.addVariable("TReturnRadius5",2.9);
  Control.addVariable("TReturnRadius6",3.5);
  Control.addVariable("TReturnRadius7",3.7);

  Control.addVariable("TReturnMat0","ParaH2");
  Control.addVariable("TReturnMat1","Aluminium");
  Control.addVariable("TReturnMat2","Void");
  Control.addVariable("TReturnMat3","Aluminium");
  Control.addVariable("TReturnMat4","Void");
  Control.addVariable("TReturnMat5","Aluminium");
  Control.addVariable("TReturnMat6","Void");
  Control.addVariable("TReturnMat7","Aluminium");

  Control.addVariable("TReturnTemp0",25.0);
  Control.addVariable("TReturnTemp1",25.0);

  Control.addVariable("TReturnActive0",3);
  Control.addVariable("TReturnActive1",15);
  Control.addVariable("TReturnActive2",63);
  Control.addVariable("TReturnActive3",255);
  //  Control.addVariable("TReturnActive4",127);


  /*  Control.addVariable("TReturnNSegIn",1);
  Control.addVariable("TReturnPPt0",Geometry::Vec3D(0, 0, 0));
  Control.addVariable("TReturnPPt1",Geometry::Vec3D(0, 30, 0));


  Control.addVariable("TReturnInRadius",1.5);
  Control.addVariable("TReturnInAlRadius",1.7);
  Control.addVariable("TReturnMidAlRadius",1.8);
  Control.addVariable("TReturnVoidRadius",2.3);
  Control.addVariable("TReturnOutAlRadius",2.5);

  Control.addVariable("TReturnInMat",1001);
  Control.addVariable("TReturnInAlMat", 13060);
  Control.addVariable("TReturnMidAlMat", 13060);
  Control.addVariable("TReturnVoidMat",0);
  Control.addVariable("TReturnOutAlMat", 13060);
  */
  // Low moderator
  Control.addVariable("LowModXStep", 0);  
  Control.addVariable("LowModYStep", 0);  
  Control.addVariable("LowModZStep",-ModZStep);
  Control.addVariable("LowModXYangle",125.15); 
  Control.addVariable("LowModZangle",0.0);
  Control.addVariable("LowModRadius",8.0);
  Control.addVariable("LowModHeight", TopModHeight);
  Control.addVariable("LowModMat",1001);
  Control.addVariable("LowModTemp",20.0);
  Control.addVariable("LowModNLayers",7);
  // al layer
  Control.addVariable("LowModHGap1",0.3);
  Control.addVariable("LowModRadGap1",0.3);
  Control.addVariable("LowModMaterial1", 13061);  // Al materk
  Control.addVariable("LowModTemp1",20.0);  
  // Vac gap
  Control.addVariable("LowModHGap2",0.5);
  Control.addVariable("LowModRadGap2",0.5);
  Control.addVariable("LowModMaterial2",0); 
  // Next Al layer
  Control.addVariable("LowModHGap3",0.25);
  Control.addVariable("LowModRadGap3",0.25);
  Control.addVariable("LowModMaterial3",13060); 
  Control.addVariable("LowModTemp3",300);//77.0);   kbat
  // He Layer
  Control.addVariable("LowModHGap4",0.001);
  Control.addVariable("LowModRadGap4",0.001);
  Control.addVariable("LowModMaterial4",0); 
  // Outer Layer
  Control.addVariable("LowModHGap5",0.25);
  Control.addVariable("LowModRadGap5",0.25);
  Control.addVariable("LowModMaterial5", 13060); 
  Control.addVariable("LowModTemp5",300.0); 
  // Clearance
  Control.addVariable("LowModHGap6",0.01);
  Control.addVariable("LowModRadGap6",0.01);
  Control.addVariable("LowModMaterial6",0); 

  Control.addVariable("LowAFlightXStep", -6.7);      // Step from centre -2.2
  Control.addVariable("LowAFlightZStep",0.0);      // Step from centre
  Control.addVariable("LowAFlightAngleXY1",30.0);  // Angle out
  Control.addVariable("LowAFlightAngleXY2",23.0);  // Angle out
  Control.addVariable("LowAFlightAngleZTop",0.0);  // Step down angle
  Control.addVariable("LowAFlightAngleZBase",0.0); // Step up angle
  Control.addVariable("LowAFlightHeight", TopModHeight);     // Full height
  Control.addVariable("LowAFlightWidth", TopFlightWidth);     // Full width 20
  Control.addVariable("LowAFlightNLiner",1);      // Liner
  Control.addVariable("LowAFlightLinerThick1",0.01);      // Liner
  Control.addVariable("LowAFlightLinerHeight1",0.01);      // Liner
  Control.addVariable("LowAFlightLinerMat1", 0);//13060);      // Liner

  Control.addVariable("LowBFlightXStep", 5.6);     // Angle 1.5
  Control.addVariable("LowBFlightZStep",0.0);      // Step from centre
  Control.addVariable("LowBFlightMasterXY",0.0);  // kbat !!! what is this ???
  Control.addVariable("LowBFlightAngleXY1",30.0);  // Angle out
  Control.addVariable("LowBFlightAngleXY2",30.0);  // Angle out
  Control.addVariable("LowBFlightAngleZTop",0.0);  // Step down angle
  Control.addVariable("LowBFlightAngleZBase",0.0); // Step up angle
  Control.addVariable("LowBFlightHeight", TopModHeight);     // Full height
  Control.addVariable("LowBFlightWidth", TopFlightWidth);     // Full width 22
  Control.addVariable("LowBFlightNLiner",1);      // Liner
  Control.addVariable("LowBFlightLinerThick1",0.01);   
  Control.addVariable("LowBFlightLinerHeight1",0.01);      // Liner
  Control.addVariable("LowBFlightLinerMat1", 0);//13060);      

  // 
  Control.addVariable("LowPreNLayers",3);   // it has to have 4 layers - otherwise problem with interseption with LowPreABlock (region 18)

  Control.addVariable("LowPreXStep", 0);
  Control.addVariable("LowPreYStep", 0);
  Control.addVariable("LowPreZStep", 0);
  Control.addVariable("LowPreXYangle", 0);
  Control.addVariable("LowPreZangle", 180); // 180

  Control.addVariable("LowPreHeight1",0.01);  // inner and butt Al shroud
  Control.addVariable("LowPreDepth1",0.01);  
  Control.addVariable("LowPreThick1",0.01);  
  Control.addVariable("LowPreMaterial1", 13060);  

  Control.addVariable("LowPreHeight2", 0.7);  // water layer below low cold mod
  Control.addVariable("LowPreDepth2",  1.7);  // water layer between target and low cold mod
  Control.addVariable("LowPreThick2",  0.7);  
  Control.addVariable("LowPreMaterial2", 1015);

  Control.addVariable("LowPreHeight3",0.2);  // outer Al shroud
  Control.addVariable("LowPreDepth3",0.2);  
  Control.addVariable("LowPreThick3",0.2);  
  Control.addVariable("LowPreMaterial3",13060);  

  Control.addVariable("LowPreHeight4",0.001);  
  Control.addVariable("LowPreDepth4",0.001);  
  Control.addVariable("LowPreThick4",0.001);  
  Control.addVariable("LowPreMaterial4",0); 
 
  Control.addVariable("LowPreNView",2);  
  Control.addVariable("LowPreViewHeight1",TopModHeight);  
  Control.addVariable("LowPreViewWidth1",14.1);  
  Control.addVariable("LowPreViewAngle1",5.0);
  Control.addVariable("LowPreViewOpenAngle1",30.0);  

  Control.addVariable("LowPreViewHeight2",TopModHeight);  
  Control.addVariable("LowPreViewWidth2",14.1);  
  Control.addVariable("LowPreViewAngle2",180.0);  
  Control.addVariable("LowPreViewOpenAngle2",30.0);  

  Control.addVariable("LowPreABlockActive",1);  
  Control.addVariable("LowPreABlockSide",0);  
  Control.addVariable("LowPreABlockWidth",4.0);  
  Control.addVariable("LowPreABlockHeight", LowPreABBlockHeight);  
  Control.addVariable("LowPreABlockLength", PreABBlockLength);  
  Control.addVariable("LowPreABlockWallThick",0.2);  
  Control.addVariable("LowPreABlockGap",0.01);  
  Control.addVariable("LowPreABlockWallMat",13060);  
  Control.addVariable("LowPreABlockWaterMat",1015);  

  // Other block
  Control.addVariable("LowPreBBlockActive",1);  
  Control.addVariable("LowPreBBlockSide",0);  
  Control.addVariable("LowPreBBlockWidth",4.0);  
  Control.addVariable("LowPreBBlockHeight", LowPreABBlockHeight);  
  Control.addVariable("LowPreBBlockLength", PreABBlockLength);
  Control.addVariable("LowPreBBlockWallThick",0.2);  
  Control.addVariable("LowPreBBlockGap",0.01);  
  Control.addVariable("LowPreBBlockWallMat", 13060);  
  Control.addVariable("LowPreBBlockWaterMat",1015);  

  // TOP MODERATOR PRE:
  Control.addVariable("TopPreNLayers", 3);

  Control.addVariable("TopPreXStep", 0);
  Control.addVariable("TopPreYStep", 0);
  Control.addVariable("TopPreZStep", 0);
  Control.addVariable("TopPreXYangle", 0);
  Control.addVariable("TopPreZangle", 180);

  Control.addVariable("TopPreHeight1",   0.01);  // inner and side Al - made very thin since there is already enough Al from TopMod
  Control.addVariable("TopPreDepth1",    0.01);  
  Control.addVariable("TopPreThick1",    0.01);  
  Control.addVariable("TopPreMaterial1", 13060);  

  Control.addVariable("TopPreHeight2", 0.7);  // thickness of water above cold mod
  Control.addVariable("TopPreDepth2", 1.7);   // thickness of water between target wheel and cold mod
  Control.addVariable("TopPreThick2",0.7);  
  Control.addVariable("TopPreMaterial2", 1015);

  Control.addVariable("TopPreHeight3",0.2);   // outer Al
  Control.addVariable("TopPreDepth3",0.2);  
  Control.addVariable("TopPreThick3",0.2);  
  Control.addVariable("TopPreMaterial3", 13060);  

  Control.addVariable("TopPreHeight4",0.001);
  Control.addVariable("TopPreDepth4",0.001);  
  Control.addVariable("TopPreThick4",0.001);  
  Control.addVariable("TopPreMaterial4", 0);

  // Control.addVariable("TopPreHeight5",0.0);  // outer void space
  // Control.addVariable("TopPreDepth5",0.0);  
  // Control.addVariable("TopPreThick5",0.0);  
  // Control.addVariable("TopPreMaterial5",0); 
 
  Control.addVariable("TopPreNView",2);  

  // opening near block A (bottom)
  Control.addVariable("TopPreViewHeight1", TopModHeight);  
  Control.addVariable("TopPreViewWidth1",14.1);  // width of the premoderator opening around theta = 150
  Control.addVariable("TopPreViewAngle1", 5);  // angular position of premoderator outer end-point around theta=100, its position also controled by TopPreViewWidth1
  Control.addVariable("TopPreViewOpenAngle1", 30.0);  // clockwise

  // opening near block B (upper)
  Control.addVariable("TopPreViewHeight2",TopModHeight);
  Control.addVariable("TopPreViewWidth2",14.1);  
  Control.addVariable("TopPreViewAngle2", 180.0);   // rotates full opening
  Control.addVariable("TopPreViewOpenAngle2", 30.0);  // controls angle of premoderator butt-end surface at the opening

  Control.addVariable("TopPreABlockActive", 1);
  Control.addVariable("TopPreABlockSide", 0);  // block location
  Control.addVariable("TopPreABlockWidth",4.0);  
  Control.addVariable("TopPreABlockHeight", TopPreABBlockHeight);  
  Control.addVariable("TopPreABlockLength", PreABBlockLength);  
  Control.addVariable("TopPreABlockWallThick",0.2);  
  Control.addVariable("TopPreABlockGap",0.01);  
  Control.addVariable("TopPreABlockWallMat", 13060);  
  Control.addVariable("TopPreABlockWaterMat",1015);  

  // Other block
  Control.addVariable("TopPreBBlockActive",   1);
  Control.addVariable("TopPreBBlockSide",     0);  
  Control.addVariable("TopPreBBlockWidth",   4.0);  
  Control.addVariable("TopPreBBlockHeight", TopPreABBlockHeight);  
  Control.addVariable("TopPreBBlockLength",  PreABBlockLength);   // 10.8
  Control.addVariable("TopPreBBlockWallThick",0.2);  
  Control.addVariable("TopPreBBlockGap",      0.01);  
  Control.addVariable("TopPreBBlockWallMat",13060);  
  Control.addVariable("TopPreBBlockWaterMat",1015);  
  //
  Control.addVariable("TopModXStep", 0);  
  Control.addVariable("TopModYStep", 0);  
  Control.addVariable("TopModZStep", ModZStep);
  Control.addVariable("TopModXYangle",54.850); 
  Control.addVariable("TopModZangle",180.0);
  Control.addVariable("TopModRadius", 8.0); // kbat
  Control.addVariable("TopModHeight", TopModHeight);
  Control.addVariable("TopModMat", 1001);
  Control.addVariable("TopModTemp",20.0);
  Control.addVariable("TopModNLayers",7);
  // al layer
  Control.addVariable("TopModHGap1",0.3);
  Control.addVariable("TopModRadGap1",0.3);
  Control.addVariable("TopModMaterial1", 13061);  // Al materk
  Control.addVariable("TopModTemp1",20.0);  
  // Vac gap
  Control.addVariable("TopModHGap2",0.5);
  Control.addVariable("TopModRadGap2",0.5);
  Control.addVariable("TopModMaterial2",0); 
  // Next Al layer
  Control.addVariable("TopModHGap3",0.25);
  Control.addVariable("TopModRadGap3",0.25);
  Control.addVariable("TopModMaterial3", 13060); 
  Control.addVariable("TopModTemp3", 300);//77.0);   // kbat
  // He Layer
  Control.addVariable("TopModHGap4",0.001);
  Control.addVariable("TopModRadGap4",0.001);
  Control.addVariable("TopModMaterial4",0); 
  // Outer Layer
  Control.addVariable("TopModHGap5",0.25);
  Control.addVariable("TopModRadGap5",0.25);
  Control.addVariable("TopModMaterial5", 13060); 
  Control.addVariable("TopModTemp5",300.0); 

  // Clearance
  Control.addVariable("TopModHGap6",0.01);
  Control.addVariable("TopModRadGap6",0.01);
  Control.addVariable("TopModMaterial6",0); 


  // TOP A FLIGHT
  Control.addVariable("TopAFlightXStep", -6.7);      // Step (shift) from centre along the y-axis (then why it's called XStep)?
  Control.addVariable("TopAFlightZStep",0.0);      // Step from centre
  // These 2 variables control flight line opening angle.
  // Sum of them gives total opening angle value.
  // They do not control premoderator opening angle.
  Control.addVariable("TopAFlightAngleXY1",30);  // Angle out (left bottom)
  Control.addVariable("TopAFlightAngleXY2",23);  // Angle out (left upper)

  Control.addVariable("TopAFlightAngleZTop",0.0);  // Step down angle
  Control.addVariable("TopAFlightAngleZBase",0.0); // Step up angle
  Control.addVariable("TopAFlightHeight", TopModHeight);     // Full height
  Control.addVariable("TopAFlightWidth", TopFlightWidth);     // Full width
  Control.addVariable("TopAFlightNLiner",1);      // Liner
  Control.addVariable("TopAFlightLinerThick1",0.01);      // Liner
  Control.addVariable("TopAFlightLinerHeight1",0.01);      // Liner
  Control.addVariable("TopAFlightLinerMat1", 0); //13060);      // Liner

  // B FLIGHT
  Control.addVariable("TopBFlightXStep", 5.6);     // Angle counterclockwise
  Control.addVariable("TopBFlightZStep",0.0);      // Step from centre

  // see comments for TopAFlightAngle
  Control.addVariable("TopBFlightAngleXY1",30.0);  // Angle out
  Control.addVariable("TopBFlightAngleXY2",30.0);  // Angle out

  Control.addVariable("TopBFlightAngleZTop",0.0);  // Step down angle
  Control.addVariable("TopBFlightAngleZBase",0.0); // Step up angle
  Control.addVariable("TopBFlightHeight",TopModHeight);     // Full height
  Control.addVariable("TopBFlightWidth", TopFlightWidth);     // Full width 20
  Control.addVariable("TopBFlightNLiner",1);      // Liner
  Control.addVariable("TopBFlightLinerThick1",0.01);   
  Control.addVariable("TopBFlightLinerHeight1",0.01);      // Liner
  Control.addVariable("TopBFlightLinerMat1", 0);// 13060);      


  Control.addVariable("ProtonBeamViewWidth", 11.5*2);  
  Control.addVariable("ProtonBeamViewHeight", 6.15*2);
  Control.addVariable("ProtonBeamLength1", 235);
  Control.addVariable("ProtonBeamRadius2", 11.5);

  Control.addVariable("ProtonBeamWindowActive", 1);  
  Control.addVariable("ProtonBeamWindowXStep",0.0);  
  Control.addVariable("ProtonBeamWindowYStep",-565.0);  // see TDR:214 "located about 4.5 m upstream of the centre of the monolith" 565 != 450 due to target offset
  Control.addVariable("ProtonBeamWindowZStep",0.0);
  Control.addVariable("ProtonBeamWindowXYangle",0.0); 
  Control.addVariable("ProtonBeamWindowZangle",0.0);
  Control.addVariable("ProtonBeamWindowLinacMat", 0); // material from the accelerator side
  Control.addVariable("ProtonBeamWindowMonolithMat", 2000); // material from the monolith side
  Control.addVariable("ProtonBeamWindowPipeRadius", 0.3); // outer radius
  Control.addVariable("ProtonBeamWindowPipeWallThick", 0.03);
  Control.addVariable("ProtonBeamWindowPipeWallMat", 13060); // TDR:215 Al-6061-T6
  Control.addVariable("ProtonBeamWindowPipeCoolingMat", 2001);
  Control.addVariable("ProtonBeamWindowFrameHeight", 30); // along z
  Control.addVariable("ProtonBeamWindowFrameWidth", 30); // along x
  Control.addVariable("ProtonBeamWindowFrameThick", 5); // along y
  Control.addVariable("ProtonBeamWindowFrameMat", 13060); // along y
  Control.addVariable("ProtonBeamWindowPanpipeHeight", 7.0); // TDR:214 7.0
  Control.addVariable("ProtonBeamWindowPanpipeWidth", 19.2); // TDR:214 19.2


  Control.addVariable("WheelShaftNLayers",3);
  Control.addVariable("WheelShaftHeight",435.0);
  Control.addVariable("WheelShaftRadius",32.0);
  Control.addVariable("WheelShaftCoolThick",1.0);
  Control.addVariable("WheelShaftCladThick",0.5);
  Control.addVariable("WheelShaftVoidThick",0.8);

  Control.addVariable("WheelCladShaftMat",26317);
  Control.addVariable("WheelMainShaftMat",26317);

  Control.addVariable("WheelXStep",0.0);  
  Control.addVariable("WheelYStep",110.0);  
  Control.addVariable("WheelZStep",0.0);
  Control.addVariable("WheelXYangle",0.0); 
  Control.addVariable("WheelZangle",0.0);
  Control.addVariable("WheelTargetHeight",8.0);
  Control.addVariable("WheelCoolantThickIn",0.15);
  Control.addVariable("WheelCoolantThickOut",0.13);
  Control.addVariable("WheelCaseThick",0.5);
  Control.addVariable("WheelVoidThick",1.0);

  Control.addVariable("WheelInnerRadius",84.27); // kbat 82
  Control.addVariable("WheelCoolantRadiusIn",114.5);
  Control.addVariable("WheelCoolantRadiusOut",124.8);
  Control.addVariable("WheelCaseRadius",125.0);
  Control.addVariable("WheelVoidRadius",126.0);

  Control.addVariable("WheelWMat", "M74001");
  Control.addVariable("WheelSteelMat",26317);
  Control.addVariable("WheelHeMat", 2002);
  Control.addVariable("WheelInnerMat",26317);

  Control.addVariable("WheelNLayers",24);
  Control.addVariable("WheelTemp", 600);

  Control.addVariable("WheelRadius1",85.65); // kbat 83.87
  Control.addVariable("WheelMatTYPE1",2); // 2002

  Control.addVariable("WheelRadius2",97.02);
  Control.addVariable("WheelMatTYPE2",3); // 74001

  Control.addVariable("WheelRadius3",97.52);
  Control.addVariable("WheelMatTYPE3",2);

  Control.addVariable("WheelRadius4",102.37);
  Control.addVariable("WheelMatTYPE4",3);

  Control.addVariable("WheelRadius5",102.77);
  Control.addVariable("WheelMatTYPE5",2);

  Control.addVariable("WheelRadius6",106.97);
  Control.addVariable("WheelMatTYPE6",3);

  Control.addVariable("WheelRadius7",107.29);
  Control.addVariable("WheelMatTYPE7",2);

  Control.addVariable("WheelRadius8",110.49);
  Control.addVariable("WheelMatTYPE8",3);

  Control.addVariable("WheelRadius9",110.78);
  Control.addVariable("WheelMatTYPE9",2);

  Control.addVariable("WheelRadius10",112.78);
  Control.addVariable("WheelMatTYPE10",3);

  Control.addVariable("WheelRadius11",113.05);
  Control.addVariable("WheelMatTYPE11",2);

  Control.addVariable("WheelRadius12",114.65);
  Control.addVariable("WheelMatTYPE12",3); 

  Control.addVariable("WheelRadius13",114.9);
  Control.addVariable("WheelMatTYPE13",2); 

  Control.addVariable("WheelRadius14",116.5);
  Control.addVariable("WheelMatTYPE14",3); 

  Control.addVariable("WheelRadius15",116.75);//
  Control.addVariable("WheelMatTYPE15",2); 

  Control.addVariable("WheelRadius16",118.35);
  Control.addVariable("WheelMatTYPE16",3); 

  Control.addVariable("WheelRadius17",118.6);
  Control.addVariable("WheelMatTYPE17",2); 

  Control.addVariable("WheelRadius18",120);
  Control.addVariable("WheelMatTYPE18",3); 

  Control.addVariable("WheelRadius19",120.25);
  Control.addVariable("WheelMatTYPE19",2); 

  Control.addVariable("WheelRadius20",121.65);
  Control.addVariable("WheelMatTYPE20",3); 

  Control.addVariable("WheelRadius21",121.9);
  Control.addVariable("WheelMatTYPE21",2); 

  Control.addVariable("WheelRadius22",123.1);
  Control.addVariable("WheelMatTYPE22",3); 

  Control.addVariable("WheelRadius23",123.35);
  Control.addVariable("WheelMatTYPE23",2); 

  Control.addVariable("WheelRadius24",124.55);
  Control.addVariable("WheelMatTYPE24",3); 

  //  Control.addVariable("WheelRadius25",124.8);
  //  Control.addVariable("WheelMatTYPE25",2); 

  // Control.addVariable("WheelRadius26",123.1);
  // Control.addVariable("WheelMatTYPE26",3); 

  // Control.addVariable("WheelRadius27",123.35);
  // Control.addVariable("WheelMatTYPE27",2); 

  // Control.addVariable("WheelRadius28",124.55);
  // Control.addVariable("WheelMatTYPE28",3); 

  // 5 meters is 500/sqrt(2) = 353.553
  Control.addVariable("F5X",  -353.553); // -218.13 // -200
  Control.addVariable("F5Y",  353.553);  // 168.39 // 568
  Control.addVariable("F5Z",  18); // must be cold moderator centre
  Control.addVariable("F5Length",  299); // ~ 5m - 2m
  Control.addVariable("F5XB",  -7.68);
  Control.addVariable("F5YB",  -2.25);
  Control.addVariable("F5ZB",  22.0);
  Control.addVariable("F5XC",  -0.62);
  Control.addVariable("F5YC",  7.98);
  Control.addVariable("F5ZC",  22.0);
  //  Control.addVariable("F5XG",  -7.68);
  //  Control.addVariable("F5YG",  -2.25);
  Control.addVariable("F5ZG",  14.00);

  Control.addVariable("F15X",  -353.553); // -218.13 // -200
  Control.addVariable("F15Y",  353.553);  // 168.39 // 568
  Control.addVariable("F15Z",  18); // must be cold moderator centre
  Control.addVariable("F15Length",  299); // ~ 5m - 2m
  Control.addVariable("F15XB",  -7.68);
  Control.addVariable("F15YB",  -2.25);
  Control.addVariable("F15ZB",  22.0);
  Control.addVariable("F15XC",  -0.62);
  Control.addVariable("F15YC",  7.98);
  Control.addVariable("F15ZC",  22.0);
  Control.addVariable("F15ZG",  14.00);

  Control.addVariable("F25X",  -353.553);
  Control.addVariable("F25Y",  353.553); 
  Control.addVariable("F25Z",  18);
  Control.addVariable("F25Length",  299);
  Control.addVariable("F25XB",  -7.68);
  Control.addVariable("F25YB",  -2.25);
  Control.addVariable("F25ZB",  22.0);
  Control.addVariable("F25XC",  -0.62);
  Control.addVariable("F25YC",  7.98);
  Control.addVariable("F25ZC",  22.0);
  Control.addVariable("F25ZG",  14.00);

  Control.addVariable("F35X",  -353.553);
  Control.addVariable("F35Y",  353.553); 
  Control.addVariable("F35Z",  18);
  Control.addVariable("F35Length",  299);
  Control.addVariable("F35XB",  -7.68);
  Control.addVariable("F35YB",  -2.25);
  Control.addVariable("F35ZB",  22.0);
  Control.addVariable("F35XC",  -0.62);
  Control.addVariable("F35YC",  7.98);
  Control.addVariable("F35ZC",  22.0);
  Control.addVariable("F35ZG",  14.00);

  Control.addVariable("TubeModXStep",0.0);  
  Control.addVariable("TubeModYStep",0.0);  
  Control.addVariable("TubeModZStep",-25.0);  
  Control.addVariable("TubeModXYangle",0.0); 
  Control.addVariable("TubeModZangle",0.0);
  Control.addVariable("TubeModWidth", 25.0);
  Control.addVariable("TubeModHeight", 20.6);
  Control.addVariable("TubeModDepth", 25.0);
  Control.addVariable("TubeModModMat", 12);
  Control.addVariable("TubeModWallMat", 13001);
  Control.addVariable("TubeModWallThick", 0.3);

  Control.addVariable("TubeModPreWidth", 2.2);
  Control.addVariable("TubeModPreDepth", 1.0);
  Control.addVariable("TubeModPreHeight", 2.2);
  Control.addVariable("TubeModPreLength", 20.0);
  Control.addVariable("TubeModPreMat", 1011);



  Control.addVariable("TopBeRefXStep", -BeRefXStep);  
  Control.addVariable("TopBeRefYStep", -BeRefYStep);  
  Control.addVariable("TopBeRefZStep", 0.0);
  Control.addVariable("TopBeRefXYangle",0.0); 
  Control.addVariable("TopBeRefZangle",0.0);
  Control.addVariable("TopBeRefRadius",29.5);
  Control.addVariable("TopBeRefHeight", 45.0-0.5); // height of Be without wall - corresponds to Alan's 45 (with wall)
  Control.addVariable("TopBeRefWallThick",0.5); // in addition to Be
  Control.addVariable("TopBeRefRefMat", "M04000");
  Control.addVariable("TopBeRefWallMat",13060);
  Control.addVariable("TopBeRefVoidThick", 0.0);

  Control.addVariable("LowBeRefXStep", BeRefXStep);  
  Control.addVariable("LowBeRefYStep", -BeRefYStep);  
  Control.addVariable("LowBeRefZStep", 0.0);
  Control.addVariable("LowBeRefXYangle",0.0); 
  Control.addVariable("LowBeRefZangle",180);
  Control.addVariable("LowBeRefRadius",29.5);
  Control.addVariable("LowBeRefHeight",45.0-0.5); // height of Be without wall - corresponds to Alan's 45 (with wall)
  Control.addVariable("LowBeRefWallThick",0.5); // in addition to Be
  Control.addVariable("LowBeRefRefMat", "M04000"); 
  Control.addVariable("LowBeRefWallMat",13060);
  Control.addVariable("LowBeRefVoidThick", 0.0);

  Control.addVariable("BulkXStep",0.0);
  Control.addVariable("BulkYStep",0.0);
  Control.addVariable("BulkZStep",0.0);
  Control.addVariable("BulkXYangle",0.0);
  Control.addVariable("BulkZangle",0.0);
  Control.addVariable("BulkNLayer",2);


  // a layer of 0 thickness (in the original geometry of SA it was clearance between Be and steel
  // since in the TDR geometry top and bottom Be reflectors are shifted with respect to each other, this clearance in any case needs to be upgraded, so I just remove it now.
  //  Control.Parse("TopBeRefRadius+TopBeRefWallThick"); // add .. cm to the distance between Be Ref and Bulk
  //  Control.addVariable("BulkRadius1");
  //  Control.addVariable("BulkHeight1",46.0);
  //   Control.addVariable("BulkDepth1",46.0);
  //   Control.addVariable("BulkMat1", 26316); // material of cell between Be and outer refle
  

  Control.addVariable("BulkRadius1",65.0);
  Control.addVariable("BulkHeight1",75.0);
  Control.addVariable("BulkDepth1",75.0);
  Control.addVariable("BulkMat1", "M26316");           // stainless

  // Bulk steel layer [No individual guides]
  Control.addVariable("BulkRadius2",200.0);
  Control.addVariable("BulkHeight2",200.0);
  Control.addVariable("BulkDepth2",200.0);
  Control.addVariable("BulkMat2", "M26005");


  // BULK FLIGHT VOID
  Control.addVariable("BulkLAFlightSideIndex",-2);   // Index
  Control.addVariable("BulkLAFlightXStep",0.0);      // Step from centre
  Control.addVariable("BulkLAFlightZStep",0.0);      // Step from centre
  Control.addVariable("BulkLAFlightAngleXY1",30.0);  // Angle out
  Control.addVariable("BulkLAFlightAngleXY2",30.0);  // Angle out
  Control.addVariable("BulkLAFlightAngleZTop",0.0);  // Step down angle
  Control.addVariable("BulkLAFlightAngleZBase",0.0); // Step up angle
  Control.addVariable("BulkLAFlightHeight",10.0);    // Full height
  Control.addVariable("BulkLAFlightWidth",23.3);     // Full width
  Control.addVariable("BulkLAFlightNLiner",0);       // Liner


  // SHUTTER BAY
  Control.addVariable("ShutterBayXStep",0.0);  
  Control.addVariable("ShutterBayYStep",0.0);  
  Control.addVariable("ShutterBayZStep",0.0);
  Control.addVariable("ShutterBayXYangle",0.0); 
  Control.addVariable("ShutterBayZangle",0.0);
  Control.addVariable("ShutterBayRadius",600.0);
  Control.addVariable("ShutterBayHeight",400.0);
  Control.addVariable("ShutterBayDepth",400.0);
  Control.addVariable("ShutterBayMat", 0); //26005); // kbat 26005);

  // Guide BAY [ All 4 same ]
  Control.addVariable("GuideBayXStep",0.0);  
  Control.addVariable("GuideBayYStep",0.0);  
  Control.addVariable("GuideBayZStep",0.0);
  Control.addVariable("GuideBayZangle",0.0);
  Control.addVariable("GuideBayViewAngle",68.0); 
  Control.addVariable("GuideBayInnerHeight",11.0);
  Control.addVariable("GuideBayInnerDepth",11.0);
  Control.addVariable("GuideBayMidRadius",170.0);
  Control.addVariable("GuideBayHeight",32.0);
  Control.addVariable("GuideBayDepth",40.0);
  Control.addVariable("GuideBayMat",26005);
  Control.addVariable("GuideBay1XYangle",0.0); 
  Control.addVariable("GuideBay2XYangle",180.0); 
  Control.addVariable("GuideBay3XYangle",0.0); 
  Control.addVariable("GuideBay4XYangle",180.0); 
  Control.addVariable("GuideBay1NItems",12);  
  Control.addVariable("GuideBay2NItems",12);  
  Control.addVariable("GuideBay3NItems",12);  
  Control.addVariable("GuideBay4NItems",12);  

  EssBeamLinesVariables(Control);
  EssConicModerator(Control);
  EssWheel(Control);


  Control.addVariable("sdefEnergy", 2000.0);  
  Control.addVariable("sdefYPos", -500.0);
  //  Control.addVariable("prdmp", "1e7 1e7");

   // STUFF FOR 90 angle:
  // Control.addVariable("LowPreABlockSide",1);  
  // Control.addVariable("LowPreBBlockSide",1);  
  // Control.addVariable("LowPreABlockActive",1);  
  // Control.addVariable("LowPreBBlockActive",1);  
  // Control.addVariable("TopPreABlockActive",1);  
  // Control.addVariable("TopPreBBlockActive",1);  

  /*
  Control.addVariable("GuideBay2XYangle",90.0); 
  Control.addVariable("GuideBay4XYangle",-90.0); 

  Control.addVariable("LowBFlightMasterXY",-90.0);  // Angle out
  Control.addVariable("LowPreViewAngle2",180.0-90.0);  
  Control.addVariable("LowPreBBlockXYangle",-90.0);  // Angle out

  Control.addVariable("LowAFlightXStep",1.5);     // Angle
  Control.addVariable("LowBFlightXStep",-1.5);     // Angle
  */
  return;
}

void
EssConicModerator(FuncDataBase& Control)
  /*!
    Create all the Conic moderator option variables
    \param Control :: DataBase
  */
{
  ELog::RegMethod RegA("essVariables[F]","EssConicModerator");
  
  // CONE MODERATOR
  Control.addVariable("LowConeModXStep",-0.0);      
  Control.addVariable("LowConeModYStep",4.0);      
  Control.addVariable("LowConeModZStep",-18.0);      
  Control.addVariable("LowConeModXYAngle",125.15);      
  Control.addVariable("LowConeModZAngle",0.0);      

  Control.addVariable("LowConeModIWidth",2.0);      
  Control.addVariable("LowConeModIHeight",2.0);      
  Control.addVariable("LowConeModOWidth",20.0);      
  Control.addVariable("LowConeModOHeight",10.0);      
  Control.addVariable("LowConeModLength",20.0);      
  Control.addVariable("LowConeModFaceThick",2.0);      
  Control.addVariable("LowConeModThick",1.5);      

  Control.addVariable("LowConeModAlThick",0.5);      

  Control.addVariable("LowConeModVacGap",0.3);      
  Control.addVariable("LowConeModWaterAlThick",0.5);      
  Control.addVariable("LowConeModWaterThick",2.0);      
  Control.addVariable("LowConeModVoidGap",0.3);      
  Control.addVariable("LowConeModWaterMat",1015);      

  Control.addVariable("LowConeModModTemp",20.0);         // Temperature of H2 
  Control.addVariable("LowConeModModMat",1001);            // Liquid H2
  Control.addVariable("LowConeModAlMat",5);              // Aluminium

  //  Control.addVariable("LowAFlightXStep",0.0);      // Step from centre


  // CONE MODERATOR
  Control.addVariable("LowConeModBXStep",-0.0);      
  Control.addVariable("LowConeModBYStep",4.0);      
  Control.addVariable("LowConeModBZStep",-18.0);      
  Control.addVariable("LowConeModBXYAngle",125.15-180.0);      
  Control.addVariable("LowConeModBZAngle",0.0);      

  Control.addVariable("LowConeModBIWidth",2.0);      
  Control.addVariable("LowConeModBIHeight",2.0);      
  Control.addVariable("LowConeModBOWidth",20.0);      
  Control.addVariable("LowConeModBOHeight",10.0);      
  Control.addVariable("LowConeModBLength",20.0);      
  Control.addVariable("LowConeModBFaceThick",2.0);      
  Control.addVariable("LowConeModBThick",1.5);      

  Control.addVariable("LowConeModBAlThick",0.5);      

  Control.addVariable("LowConeModBVacGap",0.3);      
  Control.addVariable("LowConeModBWaterAlThick",0.5);      
  Control.addVariable("LowConeModBWaterThick",2.0);      
  Control.addVariable("LowConeModBVoidGap",0.3);      
  Control.addVariable("LowConeModBWaterMat",1015);      

  Control.addVariable("LowConeModBModTemp",20.0);         // Temperature of H2 
  Control.addVariable("LowConeModBModMat",1001);            // Liquid H2
  Control.addVariable("LowConeModBAlMat", 13060);              // Aluminium

}


void
EssBeamLinesVariables(FuncDataBase& Control)
  /*!
    Create all the beamline variabes
    \param Control :: DataBase
  */
{
  ELog::RegMethod RegA("essVariables[F]","EssBeamLinesVariables");
  for(int i=0;i<4;i++)
    {
      const std::string baseKey=
	StrFunc::makeString("G",i+1)+"BLine";
      // BeamLine in guide bay
      Control.addVariable(baseKey+"XStep",0.0);  
      Control.addVariable(baseKey+"YStep",0.0);  
      Control.addVariable(baseKey+"ZStep",0.0);
      Control.addVariable(baseKey+"Zangle",0.0);
      Control.addVariable(baseKey+"Mat",1015);
      Control.addVariable(baseKey+"BeamXYAngle",0.0); 
      Control.addVariable(baseKey+"BeamZAngle",0.0);
      Control.addVariable(baseKey+"BeamXStep",0.0);
      Control.addVariable(baseKey+"BeamZStep",0.0);
      Control.addVariable(baseKey+"BeamHeight",7.2);
      Control.addVariable(baseKey+"BeamWidth",7.6);
      Control.addVariable(baseKey+"NSegment",3);
      Control.addVariable(baseKey+"Width1",22.0);
      Control.addVariable(baseKey+"Height1",22.0);
      Control.addVariable(baseKey+"Width2",28.0);
      Control.addVariable(baseKey+"Height2",44.0);
      Control.addVariable(baseKey+"Width3",40.0);
      Control.addVariable(baseKey+"Height3",60.0);
      Control.addVariable(baseKey+"Length1",170.0);
      Control.addVariable(baseKey+"Length2",170.0);
      Control.addVariable(baseKey+"1XYangle",27.50); 
      Control.addVariable(baseKey+"2XYangle",22.5); 
      Control.addVariable(baseKey+"3XYangle",17.5); 
      Control.addVariable(baseKey+"4XYangle",12.5); 
      Control.addVariable(baseKey+"5XYangle",7.5); 
      Control.addVariable(baseKey+"6XYangle",2.5); 
      Control.addVariable(baseKey+"7XYangle",-2.5); 
      Control.addVariable(baseKey+"8XYangle",-7.5);
      Control.addVariable(baseKey+"9XYangle",-12.5); 
      Control.addVariable(baseKey+"10XYangle",-17.5); 
      Control.addVariable(baseKey+"11XYangle",-22.5); 
      Control.addVariable(baseKey+"12XYangle",-27.50); 
    }
  return;
}

}  // NAMESPACE setVariable
