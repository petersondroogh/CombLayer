/********************************************************************* 
  CombLayer : MCNP(X) Input builder
 
 * File:   essBuild/linacVariables.cxx
 *
 * Copyright (c) 2017 by Konstantin Batkov
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
#include <memory>

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


// Scales are used to calculate non-marked dimensions:
// Scale of the plots in SPLTDISH0001:  1024.0/70
//const double SCALE1(14.6286);
// Scale of the plots in SPLTDISH0052 (top-left corner and central part):  620.0/106
const double SCALE52tl(5.84906);
// Scale of the plots in SPLTDISH0052 (top-right corner):  267.52/92
const double SCALE52tr(2.90783);
// Scale of the plots in SPLTDISH0052 (bottom-left corner):  14.5/10
const double SCALE52bl(1.45);
// Scale of the plots in Fc_design.pdf received from LT 13 March 2017
const double SCALE3(32.0/34.0);

namespace setVariable
{

void
EssLinacVariables(FuncDataBase& Control)
  /*!
    Create all the linac variables
    \param Control :: DataBase
  */
{
  ELog::RegMethod RegA("essVariables[F]","EssLinacVariables");

  Control.addVariable("EngineeringActive",0);

  Control.addVariable("LinacEngineeringActive",1);
  const double dtl1Start(1103.6973); // email from Carl-Johan 5 Jun 2018
  Control.addVariable("LinacLengthBack",-dtl1Start);
  Control.addVariable("LinacLengthFront",51044.4); // MARS
  Control.addVariable("LinacWidthLeft",600./2.0+15.0); // K01-20---6-G01---011
  Control.addVariable("LinacWidthRight",600./2.0-15.0); // K01-20---6-G01---011
  Control.addVariable("LinacHeight",200.0); // Height+Depth from K01-20---6-G01---011; center communicated by Lali
  Control.addVariable("LinacDepth",150.0); // Height+Depth from K01-20---6-G01---011; center communicated by Lali
  Control.addVariable("LinacWallThick",50.0); // K01-20---6-G01---011
  Control.addVariable("LinacRoofThick",80.0); // K01-20---6-G01---011
  Control.addVariable("LinacFloorThick",75.0); // K01-20---6-G01---011
  Control.addVariable("LinacFloorWidthLeft",760.0/2.0+15.0); // K01-20---6-G01---011
  Control.addVariable("LinacFloorWidthRight",760.0/2.0-15.0); // K01-20---6-G01---011
  Control.addVariable("LinacWallMat","SkanskaConcrete");
 
  Control.addVariable("LinacAirMat","Air");
  Control.addVariable("LinacSteelMat","Stainless304"); // Lali says, but promised to check
  Control.addVariable("LinacAlMat","Aluminium"); // check it !!!
  Control.addVariable("LinacWaterMat","H2O");
  Control.addVariable("LinacCopperMat","Silver%Copper%0.015"); // Copper + 0.015% weight Ag since Ag is a gamma source (YJL said)
  Control.addVariable("LinacCopperMat","Copper");
  Control.addVariable("LinacGraphiteMat","Graphite");


  // Temporary shielding walls
  // Next two lines: Lali commnuicated: 15.05.2018
  Control.addVariable("LinacNTSW", 1);
  Control.addParse<double>("LinacTSW0Length", "LinacWidthLeft+LinacWidthRight");
  // To be optimised
  Control.addVariable("LinacTSW0Width", 70.0);
  Control.addVariable("LinacTSW0XStep", 4760.0);
  Control.addVariable("LinacTSW0XYAngle", 0.0);
  Control.addVariable("LinacTSW0Mat", "SkanskaConcrete");
  Control.addVariable("LinacTSW0NLayers", 1); // for biasing

  Control.copyVarSet("LinacTSW0", "LinacTSW1");
  Control.addVariable("LinacTSW1XStep", 4110.0);
  Control.addVariable("LinacTSW1Width", 140.0);
  Control.copyVarSet("LinacTSW0", "LinacTSW2");
  Control.addVariable("LinacTSW2XStep", 4330.0);

  
  // Beam dump
  Control.addVariable("LinacBeamDumpActive", 0);
  
  Control.addVariable("LinacBeamDumpYStep",100); // arbitrary number to put if after last DTL

  // Exact definition of borated concrete needed [email 3 Nov 2016]
  Control.addVariable("LinacBeamDumpConcreteMat", "SkanskaConcrete%Boron%99");

  Control.addVariable("LinacBeamDumpFrontWallLength",4*5.0); // doc SPLTDISH0001
  Control.addVariable("LinacBeamDumpFrontWallHeight",50.0); // doc SPLTDISH0001
  Control.addVariable("LinacBeamDumpFrontWallDepth",50.0); // doc SPLTDISH0001
  Control.addVariable("LinacBeamDumpFrontWallWidth",40.0); // doc SPLTDISH0001
  Control.addVariable("LinacBeamDumpFrontWallHoleRad",16.5/2.0); // doc SPLTDISH0001

  Control.addVariable("LinacBeamDumpFrontInnerWallHeight",15.0); // doc SPLTDISH0001
  Control.addVariable("LinacBeamDumpFrontInnerWallDepth",15.0); // doc SPLTDISH0001
  Control.addVariable("LinacBeamDumpFrontInnerWallLength",2*5.0); // doc SPLTDISH0001
  Control.addVariable("LinacBeamDumpFrontInnerWallHoleRad",10.4/2.0); // doc SPLTDISH0001
  Control.addVariable("LinacBeamDumpBackInnerWallLength",2*5.0+3.0); // doc SPLTDISH0001
  //ELog::EM << "measure LinacBeamDumpBackInnerWallGapLength" << ELog::endCrit;
  Control.addVariable("LinacBeamDumpBackInnerWallGapLength",1.0); // measure!

  Control.addVariable("LinacBeamDumpSideInnerWallThick", 2*5.0); // doc SPLTDISH0001
  Control.addVariable("LinacBeamDumpSideWallThick", 10.0); // Measure SPLTDISH0001 is 0.7*SCALE1=10.24 cm, but Lali & I decided to round it to 10 cm.

  Control.addVariable("LinacBeamDumpFloorLength",100.0); // doc SPLTDISH0001
  Control.addVariable("LinacBeamDumpFloorDepth",5*5.0); // doc SPLTDISH0001
  
  Control.addVariable("LinacBeamDumpPlate25Length",35.0); // doc SPLTDISH0001
  Control.addVariable("LinacBeamDumpPlate25Depth",3.0); // doc SPLTDISH0001

  // guess based on SPLTDISH0001: from the drawing, it should be the same as
  // height of plate 15, but it looks like plate 15 height (5cm) is wrong if
  // I compare it with e.g. plate 16 height
  Control.addVariable("LinacBeamDumpPlate38Depth",3.0);

  Control.addVariable("LinacBeamDumpBackWallLength",20.0);  // email Lali 3 Nov 2016
  Control.addParse<double>("LinacBeamDumpBackWallDepth",
			   "LinacBeamDumpFrontInnerWallDepth+LinacBeamDumpFloorDepth+LinacBeamDumpPlate25Depth");

  Control.addVariable("LinacBeamDumpRoofThick",20.0); // email from Lali 3 Nov 2016
  Control.addVariable("LinacBeamDumpRoofOverhangLength", 20.0); // guess based on SPLTDISH0001
  Control.addVariable("LinacBeamDumpInnerRoofThick",2*5.0); // doc SPLTDISH0001
  
  //ELog::EM << "measure LinacBeamDumpVacPipeFrontInnerWallDist" << ELog::endCrit;
  Control.addVariable("LinacBeamDumpVacPipeFrontInnerWallDist",5.0); // measure !!!
  Control.addVariable("LinacBeamDumpVacPipeLength",63.97); // SPLTDISH0052
  Control.addVariable("LinacBeamDumpVacPipeRad",6.85*SCALE52bl/2.0); // measured on SPLTDISH0052
  Control.addVariable("LinacBeamDumpVacPipeSideWallThick",(7.55-6.85)*SCALE52bl/2.0);  // measured on SPLTDISH0052

  Control.addVariable("LinacBeamDumpVacPipeLidRmax",10.4*SCALE52bl/2.0);    // measured on SPLTDISH0052
  Control.addVariable("LinacBeamDumpVacPipeLid1Length",1.7*SCALE52bl); // measured on SPLTDISH0052
  Control.addVariable("LinacBeamDumpVacPipeLid2Length",0.65*SCALE52tl);  // measured on SPLTDISH0052
  Control.addParse<double>("LinacBeamDumpVacPipeBaseLength",
			   "LinacBeamDumpVacPipeLength-LinacBeamDumpVacPipeLid1Length-5.95"); // top-left drawing in SPLTDISH0052
  Control.addVariable("LinacBeamDumpVacPipeOuterConeOffset", 10.8*SCALE52tr); // masured on SPLTDISH0052, see pencil note in the top-rigth part
  Control.addVariable("LinacBeamDumpVacPipeInnerConeTop", 2.0*SCALE52tr); // masured on SPLTDISH0052, see pencil note in the top-rigth part
  Control.addVariable("LinacBeamDumpWallThick", 0.9*SCALE52bl);   // measured on SPLTDISH0052
  
  Control.addVariable("LinacBeamDumpWaterPipeRadius", 0.25*SCALE52tr/2.0); // measured on SPLTDISH0052
  Control.addVariable("LinacBeamDumpWaterPipeLength", 10.2*SCALE52tr); // measured on SPLTDISH0052
  Control.addVariable("LinacBeamDumpWaterPipeOffsetX", 0.22*SCALE52tl); // measured on SPLTDISH0052
  Control.addVariable("LinacBeamDumpWaterPipeOffsetZ", 10.4); // measured on SPLTDISH005
  Control.addVariable("LinacBeamDumpWaterPipeDist", 0.05*SCALE52tr); // measured on SPLTDISH005

  
  // Faraday cup
  Control.addVariable("LinacFaradayCupActive", 1.0);
  Control.addVariable("LinacFaradayCupYStep", 206.3397); // arbitrary distance from the end of last DTL
  // Dimensions are based on email from LT 13 Mar 2017 (Fc_design.pdf)
  Control.addVariable("LinacFaradayCupEngineeringActive", 1);
  Control.addVariable("LinacFaradayCupLength", 3.2);
  Control.addVariable("LinacFaradayCupOuterRadius", 3.0); // d2/2
  Control.addVariable("LinacFaradayCupInnerRadius", 2.0); // d1/2
  Control.addVariable("LinacFaradayCupFaceLength", 0.1*SCALE3);
  Control.addVariable("LinacFaradayCupFaceRadius", 2.0*SCALE3); // approx
  Control.addVariable("LinacFaradayCupAbsorberLength", 0.7);
  Control.addVariable("LinacFaradayCupAbsorberMat", "Graphite");
  Control.addVariable("LinacFaradayCupBaseLength", 1.5); // e1
  
  Control.addVariable("LinacFaradayCupCollectorLength", 0.5);
  Control.addVariable("LinacFaradayCupCollectorMat", "Graphite");

  Control.addVariable("LinacFaradayCupWallMat", "Copper");

  Control.addVariable("LinacFaradayCupNShieldLayers", 3);
  Control.addVariable("LinacFaradayCupShieldMat1","Air");
  Control.addVariable("LinacFaradayCupShieldMat2","SS304L");
  Control.addVariable("LinacFaradayCupShieldMat3","SkanskaConcrete%Boron%99");
  Control.addVariable("LinacFaradayCupShieldWidthLeft1",10.0);
  Control.addVariable("LinacFaradayCupShieldWidthLeft2",30.0);
  Control.addVariable("LinacFaradayCupShieldWidthLeft3",110.0); // 100+10 to account for LinacFaradayCupShieldWidthLeft1
  Control.addVariable("LinacFaradayCupShieldWidthRight1",10.0);
  Control.addVariable("LinacFaradayCupShieldWidthRight2",30.0);
  Control.addVariable("LinacFaradayCupShieldWidthRight3",110.0); // 100+10 to account for LinacFaradayCupShieldWidthRight1
  Control.addVariable("LinacFaradayCupShieldHeight1",10.0);
  Control.addVariable("LinacFaradayCupShieldHeight2",30.0);
  Control.addVariable("LinacFaradayCupShieldHeight3",110.0); // 100+10 to account for LinacFaradayCupShieldHeight1
  Control.addVariable("LinacFaradayCupShieldDepth1",10.0);
  Control.addVariable("LinacFaradayCupShieldDepth2",30.0);
  Control.addVariable("LinacFaradayCupShieldDepth3",120.0);
  Control.addVariable("LinacFaradayCupShieldForwardLength1",10.0);
  Control.addVariable("LinacFaradayCupShieldForwardLength2",30.0);
  Control.addVariable("LinacFaradayCupShieldForwardLength3",110.0); // 100+10 to account for air layer thick
  Control.addVariable("LinacFaradayCupShieldBackLength",170.0);

  // DTL
  const size_t nDTL = 4;
  Control.addVariable("LinacNDTLTanks", nDTL);
  // DTL lengths are from Google drive / ESS DTL
  // 02 - Mechanical development and prototype construction.pdf
  // page 21
  // PMQ parameters are from:
  // DePrisco2015: 05 - PMQ Transverse Section and Data (Lali's google drive folder)
  Control.addVariable("LinacDTL1YStep", dtl1Start);

  Control.addVariable("LinacDTL1EngineeringActive", 0);
  //  Control.addVariable("LinacDTL1Length", 768.76);
  Control.addVariable("LinacDTL1IntertankLength", 30.49389999999994-1);
  //  Control.addVariable("LinacDTL2Length", 717.11);
  Control.addVariable("LinacDTL2IntertankLength", 16.84);
  //  Control.addVariable("LinacDTL3Length", 765.26);
  Control.addVariable("LinacDTL3IntertankLength", 21.85);
  //  Control.addVariable("LinacDTL4Length", 791.71);
  Control.addVariable("LinacDTL4IntertankLength", 24.91);
  //  Control.addVariable("LinacDTL5Length", 775.79);
  Control.addVariable("LinacDTL5IntertankLength", 100.0);
  
  Control.addVariable("LinacDTLIntertankRadius", 2.8); // MARS
  Control.addVariable("LinacDTLIntertankWallThick", 0.2); // MARS


  ELog::EM << "LinacDTL1Mat4: RB: mix of Cu with water; MARS: STST. What is correct?" << ELog::endCrit;

  std::string strtmp;
  for (size_t i=1; i<=nDTL; i++)
    {
      strtmp = std::to_string(i);
      Control.addVariable("LinacDTL"+strtmp+"PMQ1Mat5", "SS304L");  // DTL1 first half

      // these variables define radii of PMQ
      Control.addVariable("LinacDTL"+strtmp+"EngineeringActive", 0);
      Control.addVariable("LinacDTL"+strtmp+"NLayers", 7);
      Control.addVariable("LinacDTL"+strtmp+"Radius1", 1.0); // DePrisco2015, table 2
      Control.addVariable("LinacDTL"+strtmp+"Mat1", "Void");
      Control.addVariable("LinacDTL"+strtmp+"Radius2", 1.15); // DePrisco2015, table 2
      Control.addVariable("LinacDTL"+strtmp+"Mat2", "Copper"); // rbfrend2-9100
      Control.addVariable("LinacDTL"+strtmp+"Radius3", 2.9);
      Control.addVariable("LinacDTL"+strtmp+"Mat3", "SS304L");
      Control.addVariable("LinacDTL"+strtmp+"Radius4", 4.5);
      Control.addVariable("LinacDTL"+strtmp+"Mat4", "SS304L");
      Control.addVariable("LinacDTL"+strtmp+"Radius5", 25.95);  // DTL_model_picture.png - email from RB 14 Mar 2017
      Control.addVariable("LinacDTL"+strtmp+"Mat5", "Void");
      Control.addVariable("LinacDTL"+strtmp+"Radius6", 26);  // MARS
      Control.addVariable("LinacDTL"+strtmp+"Mat6", "Copper");
      Control.addVariable("LinacDTL"+strtmp+"Radius7", 31);  // MARS
      Control.addVariable("LinacDTL"+strtmp+"Mat7", "SS304L");
      Control.addVariable("LinacDTL"+strtmp+"AirMat", "Air");

      Control.addVariable("LinacDTL"+strtmp+"PMQNBars", 16);      // DePrisco2015, page 1
      Control.addVariable("LinacDTL"+strtmp+"PMQBarHeight", 1.4); // DePrisco2015, fig 2
      Control.addVariable("LinacDTL"+strtmp+"PMQBarThick",  0.4); // DePrisco2015, fig 2
      Control.addVariable("LinacDTL"+strtmp+"PMQBarMat",  "Sm2Co17"); // DePrisco2015, page 1
    }

  // PMQ lengths are read from the optics file and dumped in the
  // automaticall generated pmqVariables.cxx by the optics2var script
  EssLinacPMQVariables(Control);

  // Klystron Gallery
  Control.addParse<double>("KGLengthBack", "LinacLengthBack-100"); // -100 to avoid intersection with FEB
  Control.addVariable("KGWallThick", 20);
  Control.addParse<double>("KGLengthFront", "37957.5-KGWallThick"); // ESS-0025905, number is approx
  // KG dimensions are measured from ESS-0308575 (email from Pontus 16.05)
  Control.addVariable("KGWidthLeft", 800);
  Control.addVariable("KGWidthRight", 800);
  Control.addVariable("KGHeight", 647/2.0);
  Control.addVariable("KGDepth", 647/2.0);
  Control.addVariable("KGRoofThick", 27.6);
  Control.addVariable("KGFloorThick", 20);
  Control.addVariable("KGRoofAngle", 4.0); // calculated
  Control.addVariable("KGWallMat", "SkanskaConcrete");
  Control.addVariable("KGAirMat", "Air");
  // 130 and 570 are stub legs lengths
  Control.addParse<double>("KGXStep",
			   "-(LinacWidthLeft+KGWidthRight+LinacWallThick+KGWallThick+130+570)"); // check
  Control.addVariable("KGZStep", 573.5);

  Control.addVariable("LinacNStubs", 2);
  // Stub dimensions are measured from ESS-0308575 (email from Pontus 16.05)
  // Lengths are defined as in Sullivan, page 68
  Control.addVariable("Stub100WallThick", 30); // MARS
  Control.addParse<double>("Stub100Length1",
			   "340.0+LinacWidthRight+LinacWallThick-Stub100WallThick");
  Control.addParse<double>("Stub100Length2", "610.0-2*Stub100WallThick");
  Control.addVariable("Stub100Width",   180.0);
  Control.addVariable("Stub100Height",  150.0);
  Control.addVariable("Stub100MainMat", "Air");
  Control.addVariable("Stub100WallMat", "SkanskaConcrete"); // check

  Control.addVariable("Stub100YStep", 1845.0); // ESS-0025905 and MARS

  Control.copyVarSet("Stub100", "Stub110");
  Control.addVariable("Stub110YStep", 3761.5); // ESS-0025905 and MARS

  // Ducts in the stubs

  // Stub100Duct variables according to ESS-0318785.1
  // See email from Carl-Johan 3 Aug 2018
  Control.addVariable("Stub100NDucts", 2); // 8
  Control.addVariable("Stub100Duct1WallThick", 0.0);
  Control.addVariable("Stub100Duct1Width",   70.0);
  Control.addVariable("Stub100Duct1Height",  30.0);
  Control.addParse<double>("Stub100Duct1Length1", "Stub100Length1");
  Control.addParse<double>("Stub100Duct1Length2", "400+Stub100Duct1Height");
  Control.addVariable("Stub100Duct1MainMat", "StbTCABL");
  Control.addVariable("Stub100Duct1WallMat", "Void");
  Control.addParse<double>("Stub100Duct1ZStep", "(Stub100Duct1Height-Stub100Height)/2.0");
  Control.addParse<double>("Stub100Duct1YStep", "(Stub100Duct1Width-Stub100Width)/2.0");

  Control.copyVarSet("Stub100Duct1", "Stub100Duct2");
  Control.addVariable("Stub100Duct2Height",  37.0);
  Control.addParse<double>("Stub100Duct2Length1",
			   "Stub100Length1-Stub100Height+Stub100Duct2Height");
  Control.addParse<double>("Stub100Duct2Length2", "400+Stub100Duct2Height");
  Control.addParse<double>("Stub100Duct2ZStep",
			   "(Stub100Height-Stub100Duct2Height)/2.0");
  Control.addParse<double>("Stub100Duct2YStep",
			   "(Stub100Duct2Width-Stub100Width)/2.0");

  // Duct 3-8 length are approximate
  Control.addVariable("Stub100Duct3WallThick", 0.48);
  Control.addVariable("Stub100Duct3Width",   14.61);
  Control.addVariable("Stub100Duct3Height",  58.42);
  Control.addParse<double>("Stub100Duct3Length1", "Stub100Length1-Stub100Duct1Height-20");
  Control.addParse<double>("Stub100Duct3Length2", "400+Stub100Duct3Height");
  Control.addVariable("Stub100Duct3MainMat", "Air");
  Control.addVariable("Stub100Duct3WallMat", "Aluminium");
  Control.addParse<double>("Stub100Duct3YStep", "(Stub100Duct3Width-Stub100Width)/2.0+10");

  Control.copyVarSet("Stub100Duct3", "Stub100Duct4");
  Control.addParse<double>("Stub100Duct4YStep", "Stub100Duct3YStep+Stub100Duct3Width+10");

  Control.copyVarSet("Stub100Duct4", "Stub100Duct5");
  Control.addParse<double>("Stub100Duct5YStep", "Stub100Duct4YStep+Stub100Duct4Width+75");

  Control.copyVarSet("Stub100Duct5", "Stub100Duct6");
  Control.addParse<double>("Stub100Duct6YStep", "Stub100Duct6YStep+Stub100Duct5Width+10");

  Control.copyVarSet("Stub100Duct6", "Stub100Duct7");
  Control.addVariable("Stub100Duct7Height",   14.61);
  Control.addVariable("Stub100Duct7Width",  58.42);
  Control.addParse<double>("Stub100Duct7Length1", "Stub100Duct6Length1+Stub100Duct7Height+10");
  Control.addParse<double>("Stub100Duct7Length2", "Stub100Duct6Length2-45");
  Control.addParse<double>("Stub100Duct7YStep", "-(Stub100Duct7Width-Stub100Width)/2.0-Stub100Duct7WallThick-10");
  Control.addVariable("Stub100Duct7ZStep", -50);

  Control.copyVarSet("Stub100Duct7", "Stub100Duct8");
  Control.addParse<double>("Stub100Duct8Length1", "Stub100Duct6Length1-Stub100Duct6Height-10");
  Control.addVariable("Stub100Duct8ZStep", 50);


  // Stub110Duct variables according to ESS-0318786 rev 1.0
  // See email from Carl-Johan 3 Aug 2018
  Control.addVariable("Stub110NDucts", 2);
  Control.addVariable("Stub110Duct1WallThick", 0.0);
  Control.addVariable("Stub110Duct1Width",   60.0);
  Control.addVariable("Stub110Duct1Height",  35.0);
  Control.addParse<double>("Stub110Duct1Length1", "Stub110Length1");
  Control.addParse<double>("Stub110Duct1Length2", "400+Stub110Duct1Height");
  Control.addVariable("Stub110Duct1MainMat", "StbTCABL");
  Control.addVariable("Stub110Duct1WallMat", "Void");
  Control.addParse<double>("Stub110Duct1ZStep", "(Stub110Duct1Height-Stub110Height)/2.0");
  Control.addParse<double>("Stub110Duct1YStep", "(Stub110Duct1Width-Stub110Width)/2.0");

  Control.addVariable("Stub110Duct2WallThick", 0.0);
  Control.addVariable("Stub110Duct2Width",   60.0);
  Control.addVariable("Stub110Duct2Height",  35.0);
  Control.addParse<double>("Stub110Duct2Length1", "Stub110Length1");
  Control.addParse<double>("Stub110Duct2Length2", "400+Stub110Duct2Height");
  Control.addVariable("Stub110Duct2MainMat", "StbTCABL");
  Control.addVariable("Stub110Duct2WallMat", "Void");
  Control.addParse<double>("Stub110Duct2ZStep",
			   "(Stub110Duct2Height-Stub110Height)/2.0");
  Control.addParse<double>("Stub110Duct2YStep",
			   "-(Stub110Duct2Width-Stub110Width)/2.0");

  // FEB dimensions are measured from the STEP file received from
  // Carl-Johan 31.05.2018
  Control.addVariable("FEBLength",     2000.0);
  Control.addVariable("FEBWidthRight", 2667.0);
  Control.addVariable("FEBWidthLeft",   514.0);
  Control.addVariable("FEBWallThick",    50.0);
  Control.addVariable("FEBMainMat",   "Air");
  Control.addVariable("FEBWallMat",   "SkanskaConcrete");
  Control.addVariable("FEBGapMat",   "Void"); // ???
  Control.addVariable("FEBGapALength", 180.0);
  Control.addVariable("FEBGapAHeightLow", 100.0);
  Control.addVariable("FEBGapAHeightTop", 130.0);
  Control.addVariable("FEBGapADist", 90.0);
  Control.addParse<double>("FEBGapBLength",    "FEBGapALength");
  Control.addParse<double>("FEBGapBHeightLow", "FEBGapAHeightLow");
  Control.addParse<double>("FEBGapBHeightTop", "FEBGapAHeightTop");
  Control.addParse<double>("FEBGapBDist",      "FEBGapADist");
  Control.addVariable("FEBGapABDist", 635.03);

  Control.addParse<double>("FEBShieldWall1Offset", "LinacWidthLeft");
  Control.addVariable("FEBShieldWall1Thick", 100.0);
  Control.addVariable("FEBShieldWall1Length", 1678.0);

  Control.addVariable("FEBShieldWall2Offset", 1000.0);
  Control.addParse<double>("FEBShieldWall2Thick", "FEBShieldWall1Thick");
  Control.addVariable("FEBShieldWall2Length", 800.0);

  Control.addVariable("FEBShieldWall3Offset", 715.0);
  Control.addVariable("FEBShieldWall3Thick",  100.0);
  Control.addVariable("FEBShieldWall3Length", 600.0);

  Control.addVariable("FEBDropHatchLength", 1100.0);
  Control.addVariable("FEBDropHatchWidth",   377.0);
  Control.addVariable("FEBDropHatchWallThick", 70.0);

  // RFQ
  Control.addVariable("RFQYStep", 233.697); // so the distance b/w start of RFQ and DTL1 as in rbfrend2-9102
  Control.addVariable("RFQLength", 640.0);
  Control.addVariable("RFQOuterWidth", 29.0);
  Control.addVariable("RFQInnerWidth", 18.0);
  Control.addVariable("RFQWallThick", 5.6);
  Control.addVariable("RFQMainMat", "Void");
  Control.addVariable("RFQWallMat", "Copper");
  Control.addVariable("RFQVaneThick", 2.0);

  // Berm
  Control.addParse<double>("BermLengthBack", "LinacLengthBack+100+FEBLength");
  Control.addParse<double>("BermLengthFront", "LinacLengthFront+100");
  Control.addVariable("BermHeight", 550);
  Control.addVariable("BermDepth",  500);
  Control.addVariable("BermWidthLeft", 4000);
  Control.addVariable("BermWidthRight", 3500);
  Control.addVariable("BermMat", "ClayTillLera");
  Control.addVariable("BermRoofAngle", 4.76); // approx

  return;
}

}  // NAMESPACE setVariable
