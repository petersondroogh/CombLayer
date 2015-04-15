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
#include <numeric>
#include <boost/shared_ptr.hpp>
#include <boost/array.hpp>
#include <boost/multi_array.hpp>

#include "Exception.h"
#include "FileReport.h"
#include "GTKreport.h"
#include "NameStack.h"
#include "RegMethod.h"
#include "OutputLog.h"
#include "BaseVisit.h"
#include "BaseModVisit.h"
#include "support.h"
#include "stringCombine.h"
#include "MatrixBase.h"
#include "Matrix.h"
#include "Vec3D.h"
#include "Quaternion.h"
#include "localRotate.h"
#include "masterRotate.h"
#include "Surface.h"
#include "surfIndex.h"
#include "surfRegister.h"
#include "objectRegister.h"
#include "surfEqual.h"
#include "surfDivide.h"
#include "surfDIter.h"
#include "Quadratic.h"
#include "Plane.h"
#include "Cylinder.h"
#include "Line.h"
#include "Rules.h"
#include "varList.h"
#include "Code.h"
#include "FuncDataBase.h"
#include "HeadRule.h"
#include "Object.h"
#include "Qhull.h"
#include "Simulation.h"
#include "ModelSupport.h"
#include "generateSurf.h"
#include "SimProcess.h"
#include "chipDataStore.h"
#include "LinkUnit.h"
#include "FixedComp.h"
#include "LinearComp.h"
#include "ContainedComp.h"
#include "boxValues.h"
#include "boxUnit.h"
#include "BoxLine.h"
#include "ProtonVoid.h"

namespace essSystem
{

  ProtonVoid::ProtonVoid(const std::string& Key)  :
    attachSystem::ContainedComp(),
    attachSystem::FixedComp(Key,2),
    pvIndex(ModelSupport::objectRegister::Instance().cell(Key)),
    cellIndex(pvIndex+1),protonVoidCell(0)
    /*!
      Constructor BUT ALL variable are left unpopulated.
      \param Key :: Name for item in search
    */
  {}

  ProtonVoid::ProtonVoid(const ProtonVoid& A) : 
    attachSystem::ContainedComp(A),attachSystem::FixedComp(A),
    pvIndex(A.pvIndex),cellIndex(A.cellIndex),protonVoidCell(A.protonVoidCell),
    viewWidth(A.viewWidth), viewHeight(A.viewHeight)
    /*!
      Copy constructor
      \param A :: ProtonVoid to copy
    */
  {}

  ProtonVoid&
  ProtonVoid::operator=(const ProtonVoid& A)
  /*!
    Assignment operator
    \param A :: ProtonVoid to copy
    \return *this
  */
  {
    if (this!=&A)
      {
	attachSystem::ContainedComp::operator=(A);
	attachSystem::FixedComp::operator=(A);
	cellIndex=A.cellIndex;
	protonVoidCell=A.protonVoidCell;
	viewWidth=A.viewWidth;
	viewHeight=A.viewHeight;
      }
    return *this;
  }

  ProtonVoid::~ProtonVoid() 
  /*!
    Destructor
  */
  {}

  void
  ProtonVoid::populate(const Simulation& System)
  /*!
    Populate all the variables
    \param System :: Simulation to use
  */
  {
    ELog::RegMethod RegA("ProtonVoid","populate");  
    const FuncDataBase& Control=System.getDataBase();

    // Master values
    viewWidth=Control.EvalVar<double>(keyName+"ViewWidth");
    viewHeight=Control.EvalVar<double>(keyName+"ViewHeight");
    length1=Control.EvalVar<double>(keyName+"Length1");
    radius2=Control.EvalVar<double>(keyName+"Radius2");

    return;
  }
  
  void
  ProtonVoid::createUnitVector(const attachSystem::FixedComp& FC)
  /*!
    Create the unit vectors
    - Y Down the beamline
    \param FC :: FixedComp for origin and axis
  */
  {
    ELog::RegMethod RegA("ProtonVoid","createUnitVector");

    attachSystem::FixedComp::createUnitVector(FC);
    return;
  }

  void ProtonVoid::createSurfaces()
  {
    ELog::RegMethod RegA("ProtonVoid","createSurface");

    // Void hole
    //    ModelSupport::buildCylinder(SMap,pvIndex+7,Origin,Y,viewRadius);  
    ModelSupport::buildCylinder(SMap,pvIndex+1, Origin, Y, viewWidth/2.0);

    ModelSupport::buildPlane(SMap, pvIndex+2, Origin-Z*(viewHeight/2.0), Z);
    ModelSupport::buildPlane(SMap, pvIndex+3, Origin+Z*(viewHeight/2.0), Z);
    ModelSupport::buildPlane(SMap, pvIndex+4, Origin, Y); // needed to have proton void only from one side of the wheel

    ModelSupport::buildPlane(SMap, pvIndex+5, Origin-Y*(length1), Y);
    ModelSupport::buildCylinder(SMap, pvIndex+6, Origin, Y, radius2);

    return;
  }

  void
  ProtonVoid::createObjects(Simulation& System,
			    const std::string& TargetSurfBoundary,
			    const std::string& RefSurfBoundary)
  /*!
    \param System :: Simulation to create objects in
    \param TargetSurfBoundary :: boundary layer [expect to be target edge]
    \param RefSurfBoundary :: boundary layer [expect to be reflector edge]
  */
  {
    ELog::RegMethod RegA("ProtonVoid","createObjects");
  
    std::string Out;

    Out=ModelSupport::getComposite(SMap,pvIndex, " ((-1 2 -3 -4 5) : (-6 -5)) ");
    Out+=RefSurfBoundary+" "+TargetSurfBoundary;
    System.addCell(MonteCarlo::Qhull(cellIndex++,0,0.0,Out));
    protonVoidCell=cellIndex-1;
    addOuterSurf(Out);
    //    addBoundarySurf(-SMap.realSurf(pvIndex+7));     !!! todo 

    return;
  }

  void
  ProtonVoid::createLinks()
  /*!
    Creates a full attachment set [Internal]
  */
  {
    FixedComp::setConnect(0, Origin+Y*radius2, -Y);
    FixedComp::setLinkSurf(0, SMap.realSurf(pvIndex+6));

    return;
  }

  void
  ProtonVoid::createAll(Simulation& System,
			const attachSystem::FixedComp& TargetFC,
			const long int tIndex,
			const attachSystem::FixedComp& RefFC,
			const long int rIndex)
  /*!
    Global creation of the hutch
    \param System :: Simulation to add vessel to
    \param FC :: FixedComp for origin
    \param tIndex :: Target plate surface
    \param RefFC :: FixedComp for origin
    \param rIndex :: Reflector outer surf
  */
  {
    ELog::RegMethod RegA("ProtonVoid","createAll");
    populate(System);

    createUnitVector(TargetFC);
    createSurfaces();
    const std::string TSurf=(tIndex<0) ? 
      TargetFC.getLinkComplement(static_cast<size_t>(-(tIndex+1))) : 
      TargetFC.getLinkString(static_cast<size_t>(tIndex));
    
    const std::string RSurf=(rIndex<0) ?
      RefFC.getLinkComplement(static_cast<size_t>(-(rIndex+1))) :
      RefFC.getLinkString(static_cast<size_t>(rIndex));
    createObjects(System,TSurf,RSurf);
    createLinks();
    insertObjects(System);       
    //  buildChannels(System);
    //  insertObjects(System);       

    return;
  }

void ProtonVoid::addToInsertChain(attachSystem::ContainedComp& CC) const
  /*!
    Adds this object to the containedComp to be inserted.
    \param CC :: ContainedComp object to add to this
  */
{
  for(int i=pvIndex+1;i<cellIndex;i++)
    CC.addInsertCell(i);
  
  return;
}

  
}  // NAMESPACE essSystem
