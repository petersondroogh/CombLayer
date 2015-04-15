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
#include "BeRef.h"
#include "CellBeRef.h"


namespace essSystem
{

CellBeRef::CellBeRef(const std::string& Key) :
  attachSystem::ContainedComp(),attachSystem::FixedComp(Key,3),
  refIndex(ModelSupport::objectRegister::Instance().cell(Key)),
  cellIndex(refIndex+1)
  /*!
    Constructor
    \param Key :: Name of construction key
  */
{}


CellBeRef::~CellBeRef()
  /*!
    Destructor
   */
{}

void
CellBeRef::populate(const FuncDataBase& Control)
  /*!
    Populate all the variables
    \param Control :: Variable table to use
  */
{
  ELog::RegMethod RegA("CellBeRef","populate");


    // Master values
  xStep=Control.EvalVar<double>(keyName+"XStep");
  yStep=Control.EvalVar<double>(keyName+"YStep");
  zStep=Control.EvalVar<double>(keyName+"ZStep");
  xyAngle=Control.EvalVar<double>(keyName+"XYangle");
  zAngle=Control.EvalVar<double>(keyName+"Zangle");
  radius=Control.EvalVar<double>(keyName+"Radius");   
  height=Control.EvalVar<double>(keyName+"Height");   
  wallThick=Control.EvalVar<double>(keyName+"WallThick");   

  refMat=ModelSupport::EvalMat<int>(Control,keyName+"RefMat");
  //  refMatTop=ModelSupport::EvalMat<int>(Control,keyName+"RefMatTop");
  wallMat=ModelSupport::EvalMat<int>(Control,keyName+"WallMat");   
  
  return;
}

void
CellBeRef::createUnitVector(const attachSystem::FixedComp& FC)
  /*!
    Create the unit vectors
    \param FC :: Fixed Component
  */
{
  ELog::RegMethod RegA("CellBeRef","createUnitVector");
  attachSystem::FixedComp::createUnitVector(FC);
  applyShift(xStep,yStep,zStep);
  applyAngleRotate(xyAngle,zAngle);
  
  return;
}

void
CellBeRef::createSurfaces()
  /*!
    Create Surfaces for the Be
  */
{
  ELog::RegMethod RegA("CellBeRef","createSurfaces");
      
  ModelSupport::buildCylinder(SMap,refIndex+7,Origin,Z,radius);   // SMap - local surface map [2:1300]
  ModelSupport::buildCylinder(SMap,refIndex+17,Origin,Z,radius+wallThick);  
  
  ModelSupport::buildPlane(SMap,refIndex+6,Origin+Z*(height),Z);  
  ModelSupport::buildPlane(SMap,refIndex+15, Origin+Z*(wallThick),Z); // upper/lower plane of bottom Al container (other plane is the target upper/lower surface)
  ModelSupport::buildPlane(SMap,refIndex+16, Origin+Z*(height+wallThick),Z);  

  return; 
}

void
CellBeRef::addToInsertChain(attachSystem::ContainedComp& CC) const
  /*!
    Adds this object to the containedComp to be inserted.
    \param CC :: ContainedComp object to add to this
  */
{
  for(int i=refIndex+1;i<cellIndex;i++)
    CC.addInsertCell(i);
    
  return;
}


void
CellBeRef::createObjects(Simulation& System, const std::string &TargetSurfBoundary)
  /*!
    Create the vaned moderator
    \param System :: Simulation to add results
   */
{
  ELog::RegMethod RegA("CellBeRef","createObjects");

  std::string Out;
  
  // [2:1381] There are 2 types of cells: object cells (Monte Carlo objects = MC qhulls)

  Out=ModelSupport::getComposite(SMap,refIndex," -7 -6 15");
  //  Out += TargetSurfBoundary;
  System.addCell(MonteCarlo::Qhull(cellIndex++,refMat,0.0,Out));

  //  Out=ModelSupport::getComposite(SMap,refIndex," -7 5 -6 +1 ");
  //  System.addCell(MonteCarlo::Qhull(cellIndex++,refMatTop,0.0,Out));

  Out=ModelSupport::getComposite(SMap,refIndex," -17 -16 (7:6:-15)"); // " -17 5 -16 (7:-5:6)"
  Out += TargetSurfBoundary;
  System.addCell(MonteCarlo::Qhull(cellIndex++,wallMat,0.0,Out));

  Out=ModelSupport::getComposite(SMap,refIndex," -17 -16 "); // " -17 5 -16 "
  Out += TargetSurfBoundary;
  addOuterSurf(Out); // mandatory that an object had an outer surface. it is also possible to add outer union. do not make it too big / complicated.
  return; 

}

void
CellBeRef::createLinks()
  /*!
    Creates a full attachment set
    Links/directions going outwards true.
  */
{
  /*
    Few lines of magick.
    [2:1680] How something else could connect to this object, or examine this object, or find out where it is, or someting.
    Defines how an object gets moved. Every object has to have them.
   */

  // [2:1850] The idea is that the link poin and the Z are the axis you want to connect along.
  // [2:1872] How to link/connect multiple surfacess (add...)

  FixedComp::setConnect(0,Origin+Y*radius,-Y); // origin on the y-radius (i.e. somewhere on the beam) pointing inwards (because -Y). 99% of the links are pointing outwards. Actually, -Y or +Y does not matter so much, but the surface sign does matter. Stuart says it should be +Y here.
  FixedComp::setLinkSurf(0,SMap.realSurf(refIndex+17));

  //  FixedComp::setConnect(1,Origin-Z*(height/2.0+wallThick),-Z); // down in z and the link points in negative z. This is exactly how it should be.
  //  FixedComp::setLinkSurf(1,-SMap.realSurf(refIndex+15));

  FixedComp::setConnect(2,Origin+Z*(height/2.0+wallThick),Z);
  FixedComp::setLinkSurf(2,SMap.realSurf(refIndex+16));

  return;
}


  void CellBeRef::createAll(Simulation& System, const attachSystem::FixedComp& FC, const attachSystem::FixedComp& TargetFC, const long int tIndex)
{
  /*!
    Extrenal build everything
    \param System :: Simulation
    \param FC :: FixedComponent for origin
    \param TargetFC :: target object
    \param tIndex :: link to upper/lower Target plane for upper/lower CellBeRef

    In our case the reflector has to be built relative to an origin and an axes set. 
    If you take a simple fixed object, then the axes is the axes set of this fixed object and the origin is the origin of this fixed object,
    which does not mean that the reflector and the object have the same origin, that's just the way you start and then you add the next bits.
   */

  ELog::RegMethod RegA("CellBeRef","createAll");
  // the order matters:

  populate(System.getDataBase()); // populate variables
  createUnitVector(FC); // take fixed component, then apply shift and angle rotation (transformation) for this object centre
  createSurfaces();

  const std::string TSurf=(tIndex<0) ? 
    TargetFC.getLinkComplement(static_cast<size_t>(-(tIndex+1))) : 
    TargetFC.getLinkString(static_cast<size_t>(tIndex));


  createObjects(System, TSurf);
  createLinks();
  insertObjects(System);       

  return;
}

}  // NAMESPACE instrumentSystem
