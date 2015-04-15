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
#include "Grooving.h"


namespace essSystem
{

Grooving::Grooving(const std::string& Key) :
  attachSystem::ContainedComp(),attachSystem::FixedComp(Key,5),
  grIndex(ModelSupport::objectRegister::Instance().cell(Key)),
  cellIndex(grIndex+1)
  /*!
    Constructor
    \param Key :: Name of construction key
  */
{}


Grooving::~Grooving()
  /*!
    Destructor
   */
{}

void
Grooving::populate(const FuncDataBase& Control)
  /*!
    Populate all the variables
    \param Control :: Variable table to use
  */
{
  ELog::RegMethod RegA("Grooving","populate");


    // Master values
  xStep=Control.EvalDefVar<double>(keyName+"XStep", -10);
  yStep=Control.EvalDefVar<double>(keyName+"YStep", 0);
  zStep=Control.EvalDefVar<double>(keyName+"ZStep", 15);
  xyAngle=Control.EvalDefVar<double>(keyName+"XYangle", 0);
  zAngle=Control.EvalDefVar<double>(keyName+"Zangle", 90);

  length=Control.EvalVar<double>(keyName+"Length");
  height=Control.EvalDefVar<double>(keyName+"Height", 29);   
  width=Control.EvalDefVar<double>(keyName+"Width", 2);

  wallThick=Control.EvalDefVar<double>(keyName+"WallThick", 0.2);

  nGrooves=Control.EvalDefVar<int>(keyName+"NGrooves", 2);
  dy=Control.EvalDefVar<double>(keyName+"DY", 2);

  grMat=ModelSupport::EvalDefMat<int>(Control,keyName+"GroveMat", 0);
  wallMat=ModelSupport::EvalDefMat<int>(Control,keyName+"WallMat", 13000);   

  return;
}

void
Grooving::createUnitVector(const attachSystem::FixedComp& FC)
  /*!
    Create the unit vectors
    \param FC :: Fixed Component
  */
{
  ELog::RegMethod RegA("Grooving","createUnitVector");
  attachSystem::FixedComp::createUnitVector(FC);
  applyShift(xStep,yStep,zStep);
  applyAngleRotate(xyAngle,zAngle);
  
  return;
}

void
Grooving::createSurfaces()
  /*!
    Create Surfaces for the Grooves
  */
{
  ELog::RegMethod RegA("Grooving","createSurfaces");
      
  int SI(grIndex);
  double y = -(width+wallThick*2+dy)*(nGrooves-1)/2.0;
  for (int i=0; i<nGrooves; i++) {
    ModelSupport::buildPlane(SMap, SI+2, Origin+X*(length/2+wallThick), X);  
    ModelSupport::buildPlane(SMap, SI+3, Origin-Y*(width/2+wallThick-y), Y);
    ModelSupport::buildPlane(SMap, SI+4, Origin+Y*(width/2+wallThick+y), Y);
    ModelSupport::buildPlane(SMap, SI+5, Origin-Z*(height/2+wallThick), Z);
    ModelSupport::buildPlane(SMap, SI+6, Origin+Z*(height/2+wallThick), Z);
    
    ModelSupport::buildPlane(SMap, SI+12, Origin+X*(length/2), X);  
    ModelSupport::buildPlane(SMap, SI+13, Origin-Y*(width/2-y), Y);
    ModelSupport::buildPlane(SMap, SI+14, Origin+Y*(width/2+y), Y);
    ModelSupport::buildPlane(SMap, SI+15, Origin-Z*(height/2), Z);
    ModelSupport::buildPlane(SMap, SI+16, Origin+Z*(height/2), Z);
    SI += 20;
    y += dy + width + wallThick*2;
  }
  return; 
}

void
Grooving::addToInsertChain(attachSystem::ContainedComp& CC) const
  /*!
    Adds this object to the containedComp to be inserted.
    \param CC :: ContainedComp object to add to this
  */
{
  for(int i=grIndex+1;i<cellIndex;i++)
    CC.addInsertCell(i);
    
  return;
}


void
Grooving::createObjects(Simulation& System, const std::string& s1, const std::string& s2)
  /*!
    Create the vaned moderator
    \param System :: Simulation to add results
    \param s1 :: outer surface of the moderator's shroud
    \param s2 :: inner surface of the moderator's shroud (currently not used)
   */
{
  ELog::RegMethod RegA("Grooving","createObjects");

  std::string Out;
  int SI(grIndex);
  for (int i=0; i<nGrooves; i++) {
    Out = ModelSupport::getComposite(SMap, SI, "-2 3 -4 5 -6 (12:-13:14:-15:16)") + s1;
    System.addCell(MonteCarlo::Qhull(cellIndex++, wallMat, 0.0, Out));

    Out = ModelSupport::getComposite(SMap, SI, "-12 13 -14 15 -16") + s1;
    System.addCell(MonteCarlo::Qhull(cellIndex++, grMat, 0.0, Out));

    Out = ModelSupport::getComposite(SMap, SI, "-2 3 -4 5 -6") + s1;
    addOuterUnionSurf(Out);
    SI += 20;
  }
  
  return; 

}

void
Grooving::createLinks()
  /*!
    Creates a full attachment set
    Links/directions going outwards true.
  */
{


  return;
}


void Grooving::createAll(Simulation& System, const attachSystem::FixedComp& FC, const attachSystem::FixedComp& ModFC, const long int l1, const long int l2)
{
  /*!
    Extrenal build everything
    \param System :: Simulation
    \param FC :: FixedComponent for origin
    \param ModFC :: target object
    \param l1 :: link to the outer surface of the moderator's shroud
    \param l2 :: link to the inner surface of the moderator's shroud (currently not used)
   */

  ELog::RegMethod RegA("Grooving","createAll");
  // the order matters:

  populate(System.getDataBase()); // populate variables
  createUnitVector(FC); // take fixed component, then apply shift and angle rotation (transformation) for this object centre
  createSurfaces();

  const std::string s1=(l1<0) ? ModFC.getLinkComplement(static_cast<size_t>(-(l1+1))) : ModFC.getLinkString(static_cast<size_t>(l1));
  const std::string s2=(l2<0) ? ModFC.getLinkComplement(static_cast<size_t>(-(l2+1))) : ModFC.getLinkString(static_cast<size_t>(l2));


  createObjects(System, s1, s2);
  createLinks();
  insertObjects(System);       

  return;
}

}  // NAMESPACE instrumentSystem
