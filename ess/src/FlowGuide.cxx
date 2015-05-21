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
#include <boost/shared_ptr.hpp>

#include "Exception.h"
#include "FileReport.h"
#include "GTKreport.h"
#include "NameStack.h"
#include "RegMethod.h"
#include "OutputLog.h"
#include "surfRegister.h"
#include "objectRegister.h"
#include "BaseVisit.h"
#include "BaseModVisit.h"
#include "MatrixBase.h"
#include "Matrix.h"
#include "Vec3D.h"
#include "Quaternion.h"
#include "Surface.h"
#include "surfIndex.h"
#include "Quadratic.h"
#include "Rules.h"
#include "varList.h"
#include "Code.h"
#include "FuncDataBase.h"
#include "HeadRule.h"
#include "Object.h"
#include "Qhull.h"
#include "Simulation.h"
#include "ModelSupport.h"
#include "MaterialSupport.h"
#include "generateSurf.h"
#include "support.h"
#include "stringCombine.h"
#include "LinkUnit.h"
#include "FixedComp.h"
#include "ContainedComp.h"
#include "FlowGuide.h"


namespace essSystem
{

FlowGuide::FlowGuide(const std::string& Key) :
  attachSystem::ContainedComp(),attachSystem::FixedComp(Key,3),
  onionIndex(ModelSupport::objectRegister::Instance().cell(Key)),
  cellIndex(onionIndex+1)
  /*!
    Constructor
    \param Key :: Name of construction key
  */
{
}


FlowGuide::~FlowGuide()
  /*!
    Destructor
   */
{}

void
FlowGuide::populate(const FuncDataBase& Control)
  /*!
    Populate all the variables
    \param Control :: Variable table to use
  */
{
  ELog::RegMethod RegA("FlowGuide","populate");


  // Master values
  xStep=Control.EvalDefVar<double>(keyName+"XStep", 0);
  yStep=Control.EvalDefVar<double>(keyName+"YStep", 0);
  zStep=Control.EvalDefVar<double>(keyName+"ZStep", 0);
  xyAngle=Control.EvalDefVar<double>(keyName+"XYangle", 0);
  zAngle=Control.EvalDefVar<double>(keyName+"Zangle", 0);
  //  height=Control.EvalVar<double>(keyName+"Height");   
  wallThick=Control.EvalDefVar<double>(keyName+"WallThick", 0.2);
  wallMat=ModelSupport::EvalDefMat<int>(Control,keyName+"WallMat", 0);   
  wallTemp=Control.EvalDefVar<double>(keyName+"WallTemp", 0);

  Type = Control.EvalDefVar<std::string>(keyName+"Type", "Onion");

  // variables for type 'Onion'
  nRings = Control.EvalDefVar<size_t>(keyName+"NRings", 1);
  for (size_t i=1;i<=nRings; i++) {
    radius.push_back(Control.EvalDefVar<double>(StrFunc::makeString(keyName+"Radius", i), 3));
    gateWidth.push_back(Control.EvalDefVar<double>(StrFunc::makeString(keyName+"GateWidth", i), 1));
    gateLength.push_back(Control.EvalDefVar<double>(StrFunc::makeString(keyName+"GateLength", i), 1));
  }

  // variables for type 'Star'
  // The set of variables according to the drawing received from
  // Mateusz Pucilowski 20 May 2015
  ArmYShift = Control.EvalDefVar<double>(keyName+"ArmYShift", 7);
  ArmLength = Control.EvalDefVar<double>(keyName+"ArmLength", 7.2);
  ForeArmLength = Control.EvalDefVar<double>(keyName+"ForeArmLength", 3);
  ForeArmAngle = Control.EvalDefVar<double>(keyName+"ForeArmAngle", -45);
  SideXShift = Control.EvalDefVar<double>(keyName+"SideXShift", 5.7); //??
  SideYShift = Control.EvalDefVar<double>(keyName+"SideYShift", 6); //??
  SideLength = Control.EvalDefVar<double>(keyName+"SideLength", 5.5);
  SideAngle = Control.EvalDefVar<double>(keyName+"SideAngle", -45); //??

  BaseYShift = Control.EvalDefVar<double>(keyName+"BaseYShift", ArmYShift+ArmLength-2.5); // along Y
  BaseLength = Control.EvalDefVar<double>(keyName+"BaseLength", 8.5);
  BaseArmDist = Control.EvalDefVar<double>(keyName+"BaseArmDist", wallThick/2+0.1);
  
  return;
}

void
FlowGuide::createUnitVector(const attachSystem::FixedComp& FC)
  /*!
    Create the unit vectors
    \param FC :: Fixed Component
  */
{
  ELog::RegMethod RegA("FlowGuide","createUnitVector");
  attachSystem::FixedComp::createUnitVector(FC);
  applyShift(xStep,yStep,zStep);
  applyAngleRotate(xyAngle,zAngle);
  
  return;
}

void
FlowGuide::createSurfaces()
  /*!
    Create Surfaces for the Be
  */
{
  ELog::RegMethod RegA("FlowGuide","createSurfaces");
      
  //  ModelSupport::buildPlane(SMap,onionIndex+1, Origin-Z*(height/2.0), Z);
  //  ModelSupport::buildPlane(SMap,onionIndex+2, Origin+Z*(height/2.0), Z);

  int SI(onionIndex);
  if (Type == "Onion") {
    for (size_t i=0; i<=nRings; i++) {
      ModelSupport::buildCylinder(SMap, SI+3, Origin,Z,radius[i]);
      ModelSupport::buildCylinder(SMap, SI+4,Origin,Z,radius[i]+wallThick);  
      
      // gate:
      // upper door
      ModelSupport::buildPlane(SMap, SI+5, Origin+Y*(gateWidth[i]/2.0), Y); 
      ModelSupport::buildPlane(SMap, SI+6, Origin+Y*(gateWidth[i]/2.0+wallThick), Y);
      ModelSupport::buildPlane(SMap, SI+7, Origin-X*(radius[i]+gateLength[i]), X);
      ModelSupport::buildPlane(SMap, SI+8, Origin+X*(radius[i]+gateLength[i]), X);
      
      // upper door
      ModelSupport::buildPlane(SMap, SI+9, Origin-Y*(gateWidth[i]/2.0+wallThick), Y); 
      ModelSupport::buildPlane(SMap, SI+10, Origin-Y*(gateWidth[i]/2.0), Y);
      
      SI += 10;
    }
  } else if (Type == "Star") {
    // arm
    ModelSupport::buildPlane(SMap, SI+1, Origin-X*(wallThick/2), X);
    ModelSupport::buildPlane(SMap, SI+2, Origin+X*(wallThick/2), X);
    ModelSupport::buildPlane(SMap, SI+3, Origin-Y*(ArmYShift+ArmLength), Y);
    Geometry::Vec3D yDirc(Y);
    Geometry::Quaternion::calcQRotDeg(ForeArmAngle/2, Z).rotate(yDirc);
    ModelSupport::buildPlane(SMap, SI+4, Origin-X*(wallThick/2)-Y*(ArmYShift), yDirc); // tmp

    // base
    ModelSupport::buildPlane(SMap, SI+11, Origin-Y*(BaseYShift+wallThick/2), Y);
    ModelSupport::buildPlane(SMap, SI+12, Origin-Y*(BaseYShift-wallThick/2), Y);
    //       x>0
    ModelSupport::buildPlane(SMap, SI+13, Origin+X*(BaseArmDist), X);
    ModelSupport::buildPlane(SMap, SI+14, Origin+X*(BaseArmDist+BaseLength), X);
    //       x<0
    ModelSupport::buildPlane(SMap, SI+23, Origin-X*(BaseArmDist), X);
    ModelSupport::buildPlane(SMap, SI+24, Origin-X*(BaseArmDist+BaseLength), X);

    // forearm
    yDirc = X;
    Geometry::Quaternion::calcQRotDeg(ForeArmAngle, Z).rotate(yDirc);
    ModelSupport::buildPlane(SMap, SI+31, Origin-X*(wallThick/2)-Y*(ArmYShift), yDirc); // bottom
    ModelSupport::buildPlane(SMap, SI+32, Origin+X*(wallThick/2)-Y*(ArmYShift-wallThick*tan(ForeArmAngle/2*M_PI/180)), yDirc); //top

    yDirc = Y;
    Geometry::Quaternion::calcQRotDeg(ForeArmAngle, Z).rotate(yDirc);
    ModelSupport::buildPlane(SMap, SI+33, Origin+X*(ForeArmLength*sin(-ForeArmAngle*M_PI/180.0))-Y*(ArmYShift-ForeArmLength*cos(-ForeArmAngle*M_PI/180)), yDirc);

    // sides
    //       x>0
    yDirc = X;
    Geometry::Quaternion::calcQRotDeg(SideAngle, Z).rotate(yDirc);
    ModelSupport::buildPlane(SMap, SI+41, Origin+X*(SideXShift-wallThick/2)-Y*(SideYShift), yDirc);
    ModelSupport::buildPlane(SMap, SI+42, Origin+X*(SideXShift+wallThick/2)-Y*(SideYShift-wallThick*tan(SideAngle/2*M_PI/180)), yDirc);
    yDirc = Y;
    Geometry::Quaternion::calcQRotDeg(SideAngle, Z).rotate(yDirc);
    ModelSupport::buildPlane(SMap, SI+43, Origin+X*(SideXShift-SideLength*sin(-SideAngle*M_PI/180))-Y*(SideYShift+SideLength*cos(-SideAngle*M_PI/180)), yDirc);
    ModelSupport::buildPlane(SMap, SI+44, Origin+X*(SideXShift)-Y*(SideYShift), yDirc);

    //        x<0
    yDirc = X;
    Geometry::Quaternion::calcQRotDeg(-SideAngle, Z).rotate(yDirc);
    ModelSupport::buildPlane(SMap, SI+51, Origin-X*(SideXShift-wallThick/2)-Y*(SideYShift), yDirc);
    ModelSupport::buildPlane(SMap, SI+52, Origin-X*(SideXShift+wallThick/2)-Y*(SideYShift+wallThick*tan(-SideAngle/2*M_PI/180)), yDirc);
    yDirc = Y;
    Geometry::Quaternion::calcQRotDeg(-SideAngle, Z).rotate(yDirc);
    ModelSupport::buildPlane(SMap, SI+53, Origin-X*(SideXShift+SideLength*sin(SideAngle*M_PI/180))-Y*(SideYShift+SideLength*cos(SideAngle*M_PI/180)), yDirc);
    ModelSupport::buildPlane(SMap, SI+54, Origin-X*(SideXShift)-Y*(SideYShift), yDirc);

    
  } else {
    ELog::EM << "FlowGuide: wrong type " << Type << ". Only 'Onion' or 'Star' are supported" << ELog::endErr;
  }


  return; 
}


void FlowGuide::createObjects(Simulation& System)
  /*!
    Create the onion piping
    \param System :: Simulation to add results
   */
{
  ELog::RegMethod RegA("FlowGuide","createObjects");

  std::string Out;
  
  // [2:1381] There are 2 types of cells: object cells (Monte Carlo objects = MC qhulls)

  int SI(onionIndex);
  if (Type == "Onion") {
    for (size_t i=1; i<=nRings; i++) {
      // upper half-ring
      //Out=ModelSupport::getComposite(SMap, SI, onionIndex, " (1M -2M -4) (-1M:2M:3) "); // same: (1M -2M -4 3)
      Out=ModelSupport::getComposite(SMap, SI, onionIndex, " (7 -8 5 -6 3: 6 3 -4) ") + UpperSurface + " " + BottomSurface;
      System.addCell(MonteCarlo::Qhull(cellIndex++, wallMat, wallTemp, Out));
      addOuterUnionSurf(Out);

      // bottom half-ring
      Out=ModelSupport::getComposite(SMap, SI, onionIndex, " (7 -8 -10 9 3: -9 3 -4) ") + UpperSurface + " " + BottomSurface;
      System.addCell(MonteCarlo::Qhull(cellIndex++, wallMat, wallTemp, Out));
      addOuterUnionSurf(Out);


      SI += 10;
    }
  } else if (Type=="Star") {
    // arm
    Out = ModelSupport::getComposite(SMap, SI, "1 -2 3 -4 ") + UpperSurface + " " + BottomSurface;
    System.addCell(MonteCarlo::Qhull(cellIndex++, wallMat, wallTemp, Out));
    addOuterUnionSurf(Out);

    // base
    //      x>0
    Out = ModelSupport::getComposite(SMap, SI, "11 -12 13 -14 ") + UpperSurface + " " + BottomSurface;
    System.addCell(MonteCarlo::Qhull(cellIndex++, wallMat, wallTemp, Out));
    addOuterUnionSurf(Out);
    //      x<0
    Out = ModelSupport::getComposite(SMap, SI, "11 -12 -23 24 ") + UpperSurface + " " + BottomSurface;
    System.addCell(MonteCarlo::Qhull(cellIndex++, wallMat, wallTemp, Out));
    addOuterUnionSurf(Out);

    // forearm
    Out = ModelSupport::getComposite(SMap, SI, "31 -32 4 -33 ") + UpperSurface + " " + BottomSurface;
    System.addCell(MonteCarlo::Qhull(cellIndex++, wallMat, wallTemp, Out));
    addOuterUnionSurf(Out);

    // sides
    //       x>0
    Out = ModelSupport::getComposite(SMap, SI, "41 -42 43 -44 ") + UpperSurface + " " + BottomSurface;
    System.addCell(MonteCarlo::Qhull(cellIndex++, wallMat, wallTemp, Out));
    addOuterUnionSurf(Out);

    Out = ModelSupport::getComposite(SMap, SI, "-51 52 53 -54 ") + UpperSurface + " " + BottomSurface;
    System.addCell(MonteCarlo::Qhull(cellIndex++, wallMat, wallTemp, Out));
    addOuterUnionSurf(Out);

    
  } else {
    ELog::EM << "FlowGuide: wrong type " << Type << ELog::endErr;
  }

  return; 

}

void FlowGuide::addToInsertChain(attachSystem::ContainedComp& CC) const
  /*!
    Adds this object to the containedComp to be inserted.
    \param CC :: ContainedComp object to add to this
  */
{
  for(int i=onionIndex+1; i<cellIndex; i++)
    CC.addInsertCell(i);
    
  return;
}


void
FlowGuide::createLinks()
  /*!
    Creates a full attachment set
    Links/directions going outwards true.
  */
{

  return;
}


void FlowGuide::createAll(Simulation& System, const attachSystem::FixedComp& FC)
{
  /*!
    Extrenal build everything
    \param System :: Simulation
    \param FC :: FixedComponent for origin

    In our case the reflector has to be built relative to an origin and an axes set. 
    If you take a simple fixed object, then the axes is the axes set of this fixed object and the origin is the origin of this fixed object,
    which does not mean that the reflector and the object have the same origin, that's just the way you start and then you add the next bits.
   */

  ELog::RegMethod RegA("FlowGuide","createAll");
  // the order matters:

  populate(System.getDataBase()); // populate variables
  createUnitVector(FC); // take fixed component, then apply shift and angle rotation (transformation) for this object centre
  createSurfaces();
  createObjects(System);
  createLinks();
  insertObjects(System);       

  return;
}

  void FlowGuide::setBottomSurface(const attachSystem::FixedComp& FC, const long int link)
  {
    BottomSurface  = (link<0) ? FC.getLinkComplement(static_cast<size_t>(-(link+1))) : FC.getLinkString(static_cast<size_t>(link));
    std::cout << "BottomSurface: " << BottomSurface << std::endl;
  }

  void FlowGuide::setUpperSurface(const attachSystem::FixedComp& FC, const long int link)
  {
    UpperSurface  = (link<0) ? FC.getLinkComplement(static_cast<size_t>(-(link+1))) : FC.getLinkString(static_cast<size_t>(link));
    std::cout << "UpperSurface: " << UpperSurface << std::endl;
  }


}  // NAMESPACE instrumentSystem
