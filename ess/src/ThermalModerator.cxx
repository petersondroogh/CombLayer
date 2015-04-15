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
#include "LayerComp.h"
#include "CellMap.h"
#include "ModBase.h"
#include "ThermalModerator.h"


namespace essSystem
{

  ThermalModerator::ThermalModerator(const std::string& Key) :
    constructSystem::ModBase(Key, 6)
    /*!
      Constructor
      \param Key :: Name of construction key
    */
  {}

  ThermalModerator::ThermalModerator(const ThermalModerator& A) : 
    ModBase(A)
    /*!
      Copy constructor
      \param A :: ThermalModerator to copy
    */
  {}



  ThermalModerator::~ThermalModerator()
  /*!
    Destructor
  */
  {}

  void
  ThermalModerator::populate(const FuncDataBase& Control)
  /*!
    Populate all the variables
    \param Control :: Variable table to use
  */
  {
    ELog::RegMethod RegA("ThermalModerator","populate");

    // Master values
    xStep=Control.EvalVar<double>(keyName+"XStep");
    yStep=Control.EvalVar<double>(keyName+"YStep");
    zStep=Control.EvalVar<double>(keyName+"ZStep");
    xyAngle=Control.EvalDefVar<double>(keyName+"XYangle", 0.0);
    zAngle=Control.EvalDefVar<double>(keyName+"Zangle", 0.0);
    nLayers=Control.EvalVar<size_t>(keyName+"NLayers");   

    double W, D, H, T;
    int M;
    for (size_t i=0; i<nLayers; i++) {
      if (i) {
	// *2 because container from both sides
	W += 2*Control.EvalDefVar<double>(StrFunc::makeString(keyName+"XGap", i), 0.2);
	D += 2*Control.EvalDefVar<double>(StrFunc::makeString(keyName+"YGap", i), 0.2);
	H += 2*Control.EvalDefVar<double>(StrFunc::makeString(keyName+"ZGap", i), 0.2);
	M = ModelSupport::EvalDefMat<int>(Control, StrFunc::makeString(keyName+"Mat", i), 0);
	T = (!M) ? 0.0 : Control.EvalDefVar<double>(StrFunc::makeString(keyName+"Temp", i), 0.0);
      } else {
	W = Control.EvalDefVar<double>(keyName+"Width", 0.3);
	D = Control.EvalDefVar<double>(keyName+"Depth", 0.3);
	H = Control.EvalDefVar<double>(keyName+"Height", 0.3);
	M = ModelSupport::EvalDefMat<int>(Control, keyName+"ModMat", 0);
	T = Control.EvalDefVar<double>(keyName+"Temp", 0);
      }
      //      ELog::EM << i << " " << W << " " << D << " " << H << "\t" << M << " " << T << ELog::endDiag;
      Width.push_back(W);
      Depth.push_back(D);
      Height.push_back(H);
      mat.push_back(M);
      temp.push_back(T);
    }

    WaterCylinderRad = Control.EvalVar<double>(keyName+"WaterCylinderRadius");
    WaterCylinderMat=ModelSupport::EvalDefMat<int>(Control,keyName+"WaterCylinderMat", 1011);

    AFlightAngleZBase = Control.EvalVar<double>(keyName+"AFlightAngleZBase");
    AFlightAngleZTop = Control.EvalVar<double>(keyName+"AFlightAngleZTop");
    BFlightAngleZBase = Control.EvalVar<double>(keyName+"BFlightAngleZBase");
    BFlightAngleZTop = Control.EvalVar<double>(keyName+"BFlightAngleZTop");

    AFlightAngleXY1 = Control.EvalVar<double>(keyName+"AFlightAngleXY1");
    AFlightAngleXY2 = Control.EvalVar<double>(keyName+"AFlightAngleXY2");
    BFlightAngleXY1 = Control.EvalVar<double>(keyName+"BFlightAngleXY1");
    BFlightAngleXY2 = Control.EvalVar<double>(keyName+"BFlightAngleXY2");

    FlightLineWallMat = ModelSupport::EvalDefMat<int>(Control,keyName+"FlightLineWallMat", 13000);
    FlightLineWallThick = Control.EvalDefVar<double>(keyName+"FlightLineWallThick", 0.5);
  
    return;
  }

  void
  ThermalModerator::createUnitVector(const attachSystem::FixedComp& FC)
  /*!
    Create the unit vectors
    \param FC :: Fixed Component
  */
  {
    ELog::RegMethod RegA("ThermalModerator","createUnitVector");
    attachSystem::FixedComp::createUnitVector(FC);
    applyShift(xStep,yStep,zStep);
    applyAngleRotate(xyAngle,zAngle);
  
    return;
  }

  void
  ThermalModerator::createSurfaces() 
  {
    // see archetypes/ThermalModerator
    ELog::RegMethod RegA("ThermalModerator","createSurfaces");

    double fullWidth = 0.0;
    int SI(modIndex);
    for (size_t i=0; i<nLayers; i++) {
      ModelSupport::buildPlane(SMap, SI+1, Origin-X*(Width[i]/2.0), X);
      ModelSupport::buildPlane(SMap, SI+2, Origin+X*(Width[i]/2.0), X);
      ModelSupport::buildPlane(SMap, SI+3, Origin-Y*(Depth[i]/2.0), Y);
      ModelSupport::buildPlane(SMap, SI+4, Origin+Y*(Depth[i]/2.0), Y);
      ModelSupport::buildPlane(SMap, SI+5, Origin-Z*(Height[i]/2.0), Z);
      ModelSupport::buildPlane(SMap, SI+6, Origin+Z*(Height[i]/2.0), Z);
      SI += 10;
      //      if (i==0)
	fullWidth += Width[i];
	//      else
	//	fullWidth += Width[i]*2;
    }

    // Void tube
    /// x<0
    Geometry::Vec3D yDirc(Y);
    Geometry::Quaternion::calcQRotDeg(AFlightAngleXY1, Z).rotate(yDirc);
    ModelSupport::buildPlane(SMap, SI+103, Origin-Y*(Depth[nLayers-1]/2.0 - fullWidth/4.0*tan(AFlightAngleXY1*M_PI/180.0)), yDirc);
    ModelSupport::buildPlane(SMap, SI+113, Origin-Y*(Depth[nLayers-1]/2.0 - fullWidth/4.0*tan(AFlightAngleXY1*M_PI/180.0) + FlightLineWallThick), yDirc);
    yDirc = Y;
    Geometry::Quaternion::calcQRotDeg(-AFlightAngleXY2, Z).rotate(yDirc);
    ModelSupport::buildPlane(SMap, SI+104, Origin+Y*(Depth[nLayers-1]/2.0 - fullWidth/4.0*tan(AFlightAngleXY2*M_PI/180.0)), yDirc);
    ModelSupport::buildPlane(SMap, SI+114, Origin+Y*(Depth[nLayers-1]/2.0 - fullWidth/4.0*tan(AFlightAngleXY2*M_PI/180.0) + FlightLineWallThick), yDirc);

    Geometry::Vec3D zDirc(Z);
    Geometry::Quaternion::calcQRotDeg(-AFlightAngleZBase, Y).rotate(zDirc);
    ModelSupport::buildPlane(SMap, SI+105, Origin-Z*(Height[nLayers-1]/2.0 - fullWidth/4.0*tan(AFlightAngleZBase*M_PI/180.0)), zDirc);
    ModelSupport::buildPlane(SMap, SI+115, Origin-Z*(Height[nLayers-1]/2.0 - fullWidth/4.0*tan(AFlightAngleZBase*M_PI/180.0) + FlightLineWallThick), zDirc);
    zDirc = Z;
    Geometry::Quaternion::calcQRotDeg(AFlightAngleZTop, Y).rotate(zDirc);
    ModelSupport::buildPlane(SMap, SI+106, Origin+Z*(Height[nLayers-1]/2.0 - fullWidth/4.0*tan(AFlightAngleZTop*M_PI/180.0)), zDirc);
    ModelSupport::buildPlane(SMap, SI+116, Origin+Z*(Height[nLayers-1]/2.0 - fullWidth/4.0*tan(AFlightAngleZTop*M_PI/180.0) + FlightLineWallThick), zDirc);

    // x>0
    yDirc = Y;
    Geometry::Quaternion::calcQRotDeg(-BFlightAngleXY1, Z).rotate(yDirc);
    ModelSupport::buildPlane(SMap, SI+123, Origin-Y*(Depth[nLayers-1]/2.0 - fullWidth/4.0*tan(BFlightAngleXY1*M_PI/180.0)), yDirc);
    ModelSupport::buildPlane(SMap, SI+133, Origin-Y*(Depth[nLayers-1]/2.0 - fullWidth/4.0*tan(BFlightAngleXY1*M_PI/180.0) + FlightLineWallThick), yDirc);
    yDirc = Y;
    Geometry::Quaternion::calcQRotDeg(BFlightAngleXY2, Z).rotate(yDirc);
    ModelSupport::buildPlane(SMap, SI+124, Origin+Y*(Depth[nLayers-1]/2.0 - fullWidth/4.0*tan(BFlightAngleXY2*M_PI/180.0)), yDirc);
    ModelSupport::buildPlane(SMap, SI+134, Origin+Y*(Depth[nLayers-1]/2.0 - fullWidth/4.0*tan(BFlightAngleXY2*M_PI/180.0) + FlightLineWallThick), yDirc);

    zDirc = Z;
    Geometry::Quaternion::calcQRotDeg(BFlightAngleZBase, Y).rotate(zDirc);
    ModelSupport::buildPlane(SMap, SI+125, Origin-Z*(Height[nLayers-1]/2.0 - fullWidth/4.0*tan(BFlightAngleZBase*M_PI/180.0)), zDirc);
    ModelSupport::buildPlane(SMap, SI+135, Origin-Z*(Height[nLayers-1]/2.0 - fullWidth/4.0*tan(BFlightAngleZBase*M_PI/180.0) + FlightLineWallThick), zDirc);
    zDirc = Z;
    Geometry::Quaternion::calcQRotDeg(-BFlightAngleZTop, Y).rotate(zDirc);
    ModelSupport::buildPlane(SMap, SI+126, Origin+Z*(Height[nLayers-1]/2.0 - fullWidth/4.0*tan(BFlightAngleZTop*M_PI/180.0)), zDirc);
    ModelSupport::buildPlane(SMap, SI+136, Origin+Z*(Height[nLayers-1]/2.0 - fullWidth/4.0*tan(BFlightAngleZTop*M_PI/180.0) + FlightLineWallThick), zDirc);

    // Water cylinder
    ModelSupport::buildCylinder(SMap, SI+200, Origin, Z, WaterCylinderRad);

  
    return; 
  }

  void
  ThermalModerator::addToInsertChain(attachSystem::ContainedComp& CC) const
  /*!
    Adds this object to the containedComp to be inserted.
    \param CC :: ContainedComp object to add to this
  */
  {
    for(int i=modIndex+1;i<cellIndex;i++) CC.addInsertCell(i);
    
    return;
  }


  void
  ThermalModerator::createObjects(Simulation& System,
			      const attachSystem::FixedComp& Target, const long int tIndex,
			      const attachSystem::FixedComp& ShutterBay, const long int sIndex,
			      const attachSystem::FixedComp& Reflector, const long int rSide, const long int rTop)
  {
    // see archetypes/ThermalModerator

    ELog::RegMethod RegA("ThermalModerator","createObjects");

    // target lower surface
    const std::string tSurf=(tIndex<0) ? Target.getLinkComplement(static_cast<size_t>(-(tIndex+1))) : Target.getLinkString(static_cast<size_t>(tIndex));
    //    ELog::EM << "tIndex: " << tIndex << ELog::endDiag;    ELog::EM << "tSurf: " << tSurf << ELog::endDiag;

    const std::string sSurf=(sIndex<0) ? ShutterBay.getLinkComplement(static_cast<size_t>(-(sIndex+1))) : ShutterBay.getLinkString(static_cast<size_t>(sIndex));
    //    ELog::EM << "sIndex: " << sIndex << ELog::endDiag;    ELog::EM << "sSurf: " << sSurf << ELog::endDiag;

    const std::string rSideSurf=(rSide<0) ? Reflector.getLinkComplement(static_cast<size_t>(-(rSide+1))) : Reflector.getLinkString(static_cast<size_t>(rSide));
    //    ELog::EM << "rSide: " << rSide << ELog::endDiag;    ELog::EM << "rSideSurf: " << rSideSurf << ELog::endDiag;

    const std::string rTopSurf=(rTop<0) ? Reflector.getLinkComplement(static_cast<size_t>(-(rTop+1))) : Reflector.getLinkString(static_cast<size_t>(rTop));
    ELog::EM << "rTop: " << rTop << ELog::endDiag;    ELog::EM << "rTopSurf: " << rTopSurf << ELog::endDiag;


    std::string Out;
    int SI(modIndex);

    for (size_t i=0; i<nLayers; i++) {
      Out = ModelSupport::getComposite(SMap, SI, " 1 -2 3 -4 5 -6 ");

      if (i) {
	Out += ModelSupport::getComposite(SMap, SI-10, " (-1:2:-3:4:-5:6) ");
      }

      System.addCell(MonteCarlo::Qhull(cellIndex++, mat[i], temp[i], Out));
      SI += 10;
    }

    Out = ModelSupport::getComposite(SMap, SI-10, " 1 -2 3 -4 5 -6");
    addOuterUnionSurf(Out);

    // Void tubes
    Out = ModelSupport::getComposite(SMap, SI, SI-10, " 103 -104 105 -106 -1M ") + sSurf;
    System.addCell(MonteCarlo::Qhull(cellIndex++, 0, 0.0, Out));
    Out = ModelSupport::getComposite(SMap, SI, SI-10, " 113 -114 115 -116 -1M (-103:104:-105:106:1M) ") + sSurf;
    System.addCell(MonteCarlo::Qhull(cellIndex++, FlightLineWallMat, 0.0, Out));
    Out = ModelSupport::getComposite(SMap, SI, SI-10, " 113 -114 115 -116 -1M ") + sSurf;
    addOuterUnionSurf(Out);

    Out = ModelSupport::getComposite(SMap, SI, SI-10, " 123 -124 125 -126 2M ") + sSurf;
    System.addCell(MonteCarlo::Qhull(cellIndex++, 0, 0.0, Out));
    Out = ModelSupport::getComposite(SMap, SI, SI-10, " 133 -134 135 -136 2M (-123:124:-125:126)") + sSurf;
    System.addCell(MonteCarlo::Qhull(cellIndex++, FlightLineWallMat, 0.0, Out));
    Out = ModelSupport::getComposite(SMap, SI, SI-10, " 133 -134 135 -136 2M ") + sSurf;
    addOuterUnionSurf(Out);

    
    Out = ModelSupport::getComposite(SMap, SI, SI-10, " -200 -115 -135 ( -1M:2M:-3M:4M ) ") + tSurf;
    Out += " : " + ModelSupport::getComposite(SMap, SI, SI-10, " 1M -2M 3M -4M -5M  ") + tSurf;
    addOuterUnionSurf(Out);
    System.addCell(MonteCarlo::Qhull(cellIndex++, WaterCylinderMat, 0.0, Out));
    
    return; 
  }

  void
  ThermalModerator::createLinks()
  /*!
    Creates a full attachment set
    Links/directions going outwards true.
  */
  {

    FixedComp::setConnect(0, Origin-X*(Width[nLayers-1]/2.0), -X);
    FixedComp::setLinkSurf(0, SMap.realSurf(modIndex+10*((int)nLayers-1)+1));

    FixedComp::setConnect(1, Origin-X*(Width[nLayers-2]/2.0), -X);
    FixedComp::setLinkSurf(1, SMap.realSurf(modIndex+10*((int)nLayers-2)+1));

    return;
  }

  void ThermalModerator::createAll(Simulation& System, const attachSystem::FixedComp& FC)
  {
    /*!
      This dummy funcion is here because it is required by ModBase

      Extrenal build everything
      \param System :: Simulation
      \param FC :: FixedComponent for origin

      In our case the reflector has to be built relative to an origin and an axes set. 
      If you take a simple fixed object, then the axes is the axes set of this fixed object and the origin is the origin of this fixed object,
      which does not mean that the reflector and the object have the same origin, that's just the way you start and then you add the next bits.
    */

    ELog::RegMethod RegA("ThermalModerator","createAll(Simulation &, const attachSystem::FixedComp&)");
    ELog::EM << "!!! This method should not be called !!!" << ELog::endErr;

    return;
  }



  void ThermalModerator::createAll(Simulation& System, const attachSystem::FixedComp& FC, const attachSystem::FixedComp& Target, const long int tIndex, const attachSystem::FixedComp& ShutterBay, const long int sIndex, const attachSystem::FixedComp& Reflector, const long int rSide, const long int rTop)
  {
    /*!
      Extrenal build everything
      \param System :: Simulation
      \param FC :: FixedComponent for origin

      In our case the reflector has to be built relative to an origin and an axes set. 
      If you take a simple fixed object, then the axes is the axes set of this fixed object and the origin is the origin of this fixed object,
      which does not mean that the reflector and the object have the same origin, that's just the way you start and then you add the next bits.
    */

    ELog::RegMethod RegA("ThermalModerator","createAll");
    // the order matters:

    populate(System.getDataBase()); // populate variables
    createUnitVector(FC); // take fixed component, then apply shift and angle rotation (transformation) for this object centre
    createSurfaces();
    createObjects(System, Target, tIndex, ShutterBay, sIndex, Reflector, rSide, rTop);
    createLinks();
    insertObjects(System);       

    return;
  }

  // virtual functions needed due to inheritance from LayerComp:
  // !!! copy-pasted from CylModerator.cxx, so they do not work correctly !!!
  Geometry::Vec3D
  ThermalModerator::getSurfacePoint(const size_t layerIndex,
				const size_t sideIndex) const
  /*!
    Given a side and a layer calculate the link point
    \param sideIndex :: Side [0-5]
    \param layerIndex :: layer, 0 is inner moderator [0-6]
    \return Surface point
  */
  {
    ELog::RegMethod RegA("ThermalModerator","getSurfacePoint");

    if (sideIndex>5) 
      throw ColErr::IndexError<size_t>(sideIndex,5,"sideIndex ");
    if (layerIndex>nLayers) 
      throw ColErr::IndexError<size_t>(layerIndex,nLayers,"layer");

    // Modification map:
    //  std::cout << "layerIndex sideIndex: " << layerIndex << " " << sideIndex << std::endl;
    switch(sideIndex)
      {
      case 0:
	return 0; // kbat: see CylModerator.cxx
      }
    throw ColErr::IndexError<size_t>(sideIndex,5,"sideIndex ");
  }

  std::string
  ThermalModerator::getLayerString(const size_t layerIndex,
			       const size_t sideIndex) const
  /*!
    Given a side and a layer calculate the link surf
    \param sideIndex :: Side [0-5]
    \param layerIndex :: layer, 0 is inner moderator [0-4]
    \return Surface string
  */
  {
    ELog::RegMethod RegA("ThermalModerator","getLayerString");

    if (layerIndex>nLayers) 
      throw ColErr::IndexError<size_t>(layerIndex,nLayers,"layer");

    const int SI(modIndex+static_cast<int>(layerIndex)*10);
    std::ostringstream cx;
    switch(sideIndex)
      {
      case 0:
	return ModelSupport::getComposite(SMap,SI,modIndex," 7 -2M ");
      case 1:
	return ModelSupport::getComposite(SMap,SI,modIndex," 7 2M ");
      case 2:
	return ModelSupport::getComposite(SMap,SI,modIndex," 7 -1M ");
      case 3:
	return ModelSupport::getComposite(SMap,SI,modIndex," 7 1M ");
      case 4:
	cx<<" "<<-SMap.realSurf(SI+5)<<" ";
	return cx.str();
      case 5:
	cx<<" "<<SMap.realSurf(SI+6)<<" ";
	return cx.str();
      }
    throw ColErr::IndexError<size_t>(sideIndex,5,"sideIndex ");
  }

  int
  ThermalModerator::getLayerSurf(const size_t layerIndex,
			     const size_t sideIndex) const
  /*!
    Given a side and a layer calculate the link surf
    \param sideIndex :: Side [0-5]
    \param layerIndex :: layer, 0 is inner moderator [0-4]
    \return Surface string
  */
  {
    ELog::RegMethod RegA("ThermalModerator","getLayerSurf");

    if (layerIndex>nLayers) 
      throw ColErr::IndexError<size_t>(layerIndex,nLayers,"layer");
  
    const int SI(modIndex+static_cast<int>(layerIndex)*10);
    switch(sideIndex)
      {
      case 0:
	return SMap.realSurf(SI+7);
      case 1:
	return SMap.realSurf(SI+7);
      case 2:
	return SMap.realSurf(SI+7);
      case 3:
	return SMap.realSurf(SI+7);
      case 4:
	return -SMap.realSurf(SI+5);
      case 5:
	return SMap.realSurf(SI+6);
      }
    throw ColErr::IndexError<size_t>(sideIndex,5,"sideIndex ");
  }

  int
  ThermalModerator::getCommonSurf(const size_t sideIndex) const
  /*!
    Given a side calculate the boundary surface
    \param sideIndex :: Side [0-5]
    \return Common dividing surface [outward pointing]
  */
  {
    ELog::RegMethod RegA("ThermalModerator","getCommonSurf");

    switch(sideIndex)
      {
      case 0:
	return -SMap.realSurf(modIndex+2);
      case 1:
	return SMap.realSurf(modIndex+2);
      case 2:
	return -SMap.realSurf(modIndex+1);
      case 3:
	return SMap.realSurf(modIndex+1);
      case 4:
      case 5:
	return 0;
      }
    throw ColErr::IndexError<size_t>(sideIndex,5,"sideIndex ");
  }

  ThermalModerator*
  ThermalModerator::clone() const
  /*!
    Clone copy constructor
    \return copy of this
  */
  {
    return new ThermalModerator(*this);
  }


}  // NAMESPACE instrumentSystem
