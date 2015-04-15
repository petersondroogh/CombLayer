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
#include "RecModerator.h"


namespace essSystem
{

  RecModerator::RecModerator(const std::string& Key) :
    constructSystem::ModBase(Key, 6)
    /*!
      Constructor
      \param Key :: Name of construction key
    */
  {}

  RecModerator::RecModerator(const RecModerator& A) : 
    ModBase(A)
    /*!
      Copy constructor
      \param A :: RecModerator to copy
    */
  {}



  RecModerator::~RecModerator()
  /*!
    Destructor
  */
  {}

  void
  RecModerator::populate(const FuncDataBase& Control)
  /*!
    Populate all the variables
    \param Control :: Variable table to use
  */
  {
    ELog::RegMethod RegA("RecModerator","populate");

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
      ELog::EM << i << " " << W << " " << D << " " << H << "\t" << M << " " << T << ELog::endDiag;
      Width.push_back(W);
      Depth.push_back(D);
      Height.push_back(H);
      mat.push_back(M);
      temp.push_back(T);
    }

    windowThick = Control.EvalDefVar<double>(keyName+"WindowThick", (Width[nLayers-1]-Width[nLayers-2])/4.0);
    //    ELog::EM << "windowThick: " << windowThick << ELog::endDiag;

    preWidth=Control.EvalVar<double>(keyName+"PreWidth");   
    preDepth=Control.EvalVar<double>(keyName+"PreDepth");   
    preHeight1=Control.EvalVar<double>(keyName+"PreHeight1");   
    preHeight2=Control.EvalVar<double>(keyName+"PreHeight2");   

    preWingThick=Control.EvalVar<double>(keyName+"PreWingThick");   

    preMat=ModelSupport::EvalMat<int>(Control,keyName+"PreMat");
    preTemp=Control.EvalDefVar<double>(keyName+"PreTemp", 0);
    preWallMat=ModelSupport::EvalMat<int>(Control,keyName+"PreWallMat");
    preWallThick=Control.EvalVar<double>(keyName+"PreWallThick");

    tubeHeight=Control.EvalVar<double>(keyName+"TubeHeight");
    tubeDepth=Control.EvalVar<double>(keyName+"TubeDepth");

    targetGapMat=ModelSupport::EvalDefMat<int>(Control,keyName+"TargetGapMat", 0);
    fastRefMat=ModelSupport::EvalDefMat<int>(Control,keyName+"FastRefMat", 82000);

    WaterCylinderRad = Control.EvalVar<double>(keyName+"WaterCylinderRadius");
    WaterCylinderMat=ModelSupport::EvalDefMat<int>(Control,keyName+"WaterCylinderMat", 1011);
  
    return;
  }

  void
  RecModerator::createUnitVector(const attachSystem::FixedComp& FC)
  /*!
    Create the unit vectors
    \param FC :: Fixed Component
  */
  {
    ELog::RegMethod RegA("RecModerator","createUnitVector");
    attachSystem::FixedComp::createUnitVector(FC);
    applyShift(xStep,yStep,zStep);
    applyAngleRotate(xyAngle,zAngle);
  
    return;
  }

  void
  RecModerator::createSurfaces() 
  {
    // see archetypes/RecModerator
    ELog::RegMethod RegA("RecModerator","createSurfaces");

    int SI(modIndex);
    for (size_t i=0; i<nLayers; i++) {
      ModelSupport::buildPlane(SMap, SI+1, Origin-X*(Width[i]/2.0), X);
      ModelSupport::buildPlane(SMap, SI+2, Origin+X*(Width[i]/2.0), X);
      ModelSupport::buildPlane(SMap, SI+3, Origin-Y*(Depth[i]/2.0), Y);
      ModelSupport::buildPlane(SMap, SI+4, Origin+Y*(Depth[i]/2.0), Y);
      ModelSupport::buildPlane(SMap, SI+5, Origin-Z*(Height[i]/2.0), Z);
      ModelSupport::buildPlane(SMap, SI+6, Origin+Z*(Height[i]/2.0), Z);
      SI += 10;
    }

    // premoderator
    /// upper part
    //// water
    //// x-planes are from void of the layered structure...
    ModelSupport::buildPlane(SMap, SI+3, Origin-Y*(tubeDepth/2.0-preWallThick), Y);
    ModelSupport::buildPlane(SMap, SI+4, Origin+Y*(tubeDepth/2.0-preWallThick), Y);
    // bottom Z is the upper plane in the layered structure
    ModelSupport::buildPlane(SMap, SI+5, Origin-Z*(Height[nLayers-1]/2.0+preHeight1), Z);
    //// aluminum container
    //// x-planes are from Al of the layered structure...
    ModelSupport::buildPlane(SMap, SI+15, Origin-Z*(Height[nLayers-1]/2.0+preHeight1+preWallThick), Z);

   
    /// lower part
    //// horizontal water
    ModelSupport::buildPlane(SMap, SI+13, Origin-Y*(Depth[nLayers-1]/2.0+preDepth), Y);
    ModelSupport::buildPlane(SMap, SI+14, Origin+Y*(Depth[nLayers-1]/2.0+preDepth), Y);
    ModelSupport::buildPlane(SMap, SI+26, Origin+Z*(Height[nLayers-1]/2.0+preHeight2), Z);
    //// vertical water wings
    ///// water - x and y planes are the inside ones
    ModelSupport::buildPlane(SMap, SI+31, Origin-X*(Width[nLayers-1]/2.0-preWingThick), X);
    ModelSupport::buildPlane(SMap, SI+32, Origin+X*(Width[nLayers-1]/2.0-preWingThick), X);
    const double voidZmin = tubeHeight-(Height[nLayers-1]/2.0+preHeight1+preWallThick);
    ModelSupport::buildPlane(SMap, SI+36, Origin+Z*(voidZmin-preWallThick), Z);

    // neutron extraction window
    ModelSupport::buildPlane(SMap, SI+41, Origin-X*(Width[nLayers-2]/2.0+windowThick), X);
    ModelSupport::buildPlane(SMap, SI+42, Origin+X*(Width[nLayers-2]/2.0+windowThick), X);

    
    // Void tube
    //    ModelSupport::buildPlane(SMap, SI+101, Origin-X*(500), X); // temporary
    //    ModelSupport::buildPlane(SMap, SI+102, Origin+X*(500), X); // temporary
    ModelSupport::buildPlane(SMap, SI+103, Origin-Y*(tubeDepth/2.0), Y);
    ModelSupport::buildPlane(SMap, SI+104, Origin+Y*(tubeDepth/2.0), Y);
    ModelSupport::buildPlane(SMap, SI+106, Origin+Z*(voidZmin), Z);

    // cutting by inclined plane
    Geometry::Vec3D xDirc(Z);
    Geometry::Quaternion::calcQRotDeg(27.0,X).rotate(xDirc);
    ModelSupport::buildPlane(SMap,SI+111,Origin+Z*8.0, xDirc);

    // Water cylinder
    ModelSupport::buildCylinder(SMap, SI+200, Origin, Z, WaterCylinderRad);

  
    return; 
  }

  void
  RecModerator::addToInsertChain(attachSystem::ContainedComp& CC) const
  /*!
    Adds this object to the containedComp to be inserted.
    \param CC :: ContainedComp object to add to this
  */
  {
    for(int i=modIndex+1;i<cellIndex;i++) CC.addInsertCell(i);
    
    return;
  }


  void
  RecModerator::createObjects(Simulation& System,
			      const attachSystem::FixedComp& Target, const long int tIndex,
			      const attachSystem::FixedComp& ShutterBay, const long int sIndex,
			      const attachSystem::FixedComp& Reflector, const long int rSide, const long int rTop)
  {
    // see archetypes/RecModerator

    ELog::RegMethod RegA("RecModerator","createObjects");

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

      if ((i+1)==nLayers) {
	//addOuterSideSurf(Out); do not need this since in any case full moderator is inside water premoderator
	// remove neutron extraction windows
	Out += ModelSupport::getComposite(SMap, SI, modIndex," (-1:51:-3M:4M:-5M:6M) ");
	Out += ModelSupport::getComposite(SMap, SI, modIndex," (2:-52:-3M:4M:-5M:6M) ");
      }
      if (i) {
	Out += ModelSupport::getComposite(SMap, SI-10, " (-1:2:-3:4:-5:6) ");
      }

      System.addCell(MonteCarlo::Qhull(cellIndex++, mat[i], temp[i], Out));
      SI += 10;
    }

    // neutron extraction windows (regions of thinner Al container)
    Out  = ModelSupport::getComposite(SMap, SI-10, " 1 -51 "); // 51(SI-10) = 41(SI)
    Out += ModelSupport::getComposite(SMap, SI-10, modIndex, " 3M -4M 5M -6M ");
    System.addCell(MonteCarlo::Qhull(cellIndex++, 0, 0.0, Out));

    Out  = ModelSupport::getComposite(SMap, SI-10, " -2 52 "); // 52(SI-10) = 42(SI)
    Out += ModelSupport::getComposite(SMap, SI-10, modIndex, " 3M -4M 5M -6M ");
    System.addCell(MonteCarlo::Qhull(cellIndex++, 0, 0.0, Out));
  
    // premoderator
    // water around moderator
    Out  = ModelSupport::getComposite(SMap, SI-20, SI, " 1 -2 13M -14M " );
    Out += ModelSupport::getComposite(SMap, SI, " 5 -26 " );
    addOuterUnionSurf(Out);
    Out += ModelSupport::getComposite(SMap, SI-10, " (-1:2:-3:4:-5:6) "); // complement of moderator outer layer
    System.addCell(MonteCarlo::Qhull(cellIndex++, preMat, preTemp, Out));

    // water wings
    Out  = ModelSupport::getComposite(SMap, SI, SI-20, " -2M 32 3 -4 5 -36 ");
    Out += ModelSupport::getComposite(SMap, SI, SI-10, " (-13:14:-5:26) ");
    addOuterUnionSurf(Out);
    System.addCell(MonteCarlo::Qhull(cellIndex++, preMat, preTemp, Out));

    Out  = ModelSupport::getComposite(SMap, SI, SI-20, " 1M -31 3 -4 5 -36 ");
    Out += ModelSupport::getComposite(SMap, SI, SI-10, " (-13:14:-5:26) ");
    addOuterUnionSurf(Out);
    System.addCell(MonteCarlo::Qhull(cellIndex++, preMat, preTemp, Out));

    // Al structure
    Out  = ModelSupport::getComposite(SMap, SI, SI-10, " 1M -2M 103 -104 15 -106 ");
    Out  += ModelSupport::getComposite(SMap, SI, SI-20, " (-1M:2M:-3:4:-5:36) ");
    addOuterUnionSurf(Out);
    Out  += ModelSupport::getComposite(SMap, SI-10, " (-1:2:-3:4:-5:6) ");
    System.addCell(MonteCarlo::Qhull(cellIndex++, preWallMat, preTemp, Out));
    
    // Void tube
    Out = ModelSupport::getComposite(SMap, SI, SI-10, " 2M 103 -104 15 -106 ") + sSurf;
    addOuterUnionSurf(Out);
    System.addCell(MonteCarlo::Qhull(cellIndex++, 0, 0.0, Out));

    Out = ModelSupport::getComposite(SMap, SI, SI-10, " -1M 103 -104 15 -106 ") + sSurf;
    addOuterUnionSurf(Out);
    System.addCell(MonteCarlo::Qhull(cellIndex++, 0, 0.0, Out));

    /*
    // Void/steel under the target
    //    Out = ModelSupport::getComposite(SMap, SI-10, SI, " -15M -4M") + tSurf + rSideSurf;
    Out = ModelSupport::getComposite(SMap, SI-10, SI, "1 -2 13M -14M -15M -4M") + tSurf + rSideSurf;
    addOuterUnionSurf(Out);
    System.addCell(MonteCarlo::Qhull(cellIndex++, targetGapMat, 0.0, Out));

    // Part of BeRef  under the target is a fast reflector:
    /// region inside Al box and between water wings
    // uncomment it and comment next if inclined plane is not necessary
    //    Out = ModelSupport::getComposite(SMap, SI, " 31 -32 -13 3 -36 5 ");  addOuterUnionSurf(Out); System.addCell(MonteCarlo::Qhull(cellIndex++, fastRefMat, 0.0, Out));
    // divide it by inclined plane
    Out = ModelSupport::getComposite(SMap, SI, " 31 -32 -13 3 -36 111 "); addOuterUnionSurf(Out); System.addCell(MonteCarlo::Qhull(cellIndex++, fastRefMat, 0.0, Out));
    Out = ModelSupport::getComposite(SMap, SI, " 31 -32 -13 3 -111 5 "); addOuterUnionSurf(Out); System.addCell(MonteCarlo::Qhull(cellIndex++, 0, 0.0, Out));

    // region outside Al box, same z
    //    Out = ModelSupport::getComposite(SMap, SI, " -103 -36 5 ") + rSideSurf; addOuterUnionSurf(Out); System.addCell(MonteCarlo::Qhull(cellIndex++, fastRefMat, 0.0, Out));
    // divide it by inclined plane
    Out = ModelSupport::getComposite(SMap, SI, " -103 -36 111 ") + rSideSurf + tSurf; addOuterUnionSurf(Out); System.addCell(MonteCarlo::Qhull(cellIndex++, fastRefMat, 0.0, Out));
    Out = ModelSupport::getComposite(SMap, SI, " -103 -36 5 -111 ") + rSideSurf; addOuterUnionSurf(Out); System.addCell(MonteCarlo::Qhull(cellIndex++, 0, 0.0, Out));

    // region outside Al box, z up
    // without inclined plane
    //Out = ModelSupport::getComposite(SMap, SI, " -103 -5 ") + rSideSurf + tSurf; addOuterUnionSurf(Out); System.addCell(MonteCarlo::Qhull(cellIndex++, fastRefMat, 0.0, Out));
    // with inclined plane
    Out = ModelSupport::getComposite(SMap, SI, " -103 -5 -111") + rSideSurf + tSurf; addOuterUnionSurf(Out); System.addCell(MonteCarlo::Qhull(cellIndex++, 0, 0.0, Out));
    //    Out = ModelSupport::getComposite(SMap, SI, " -103 -5 111 ") + rSideSurf + tSurf; addOuterUnionSurf(Out); System.addCell(MonteCarlo::Qhull(cellIndex++, fastRefMat, 0.0, Out));


    // region above the region inside Al box
    Out = ModelSupport::getComposite(SMap, SI, " -13 3 -15 ") + tSurf + rSideSurf;
    addOuterUnionSurf(Out);
    System.addCell(MonteCarlo::Qhull(cellIndex++, 0, 0.0, Out));
    // region outside Al box, z bot
    Out = ModelSupport::getComposite(SMap, SI, " -103  36 ") + rSideSurf + rTopSurf;
    addOuterUnionSurf(Out);
    System.addCell(MonteCarlo::Qhull(cellIndex++, fastRefMat, 0.0, Out));
    // region inside Al box, z bot
    Out = ModelSupport::getComposite(SMap, SI, " -13 3 106 ") + rTopSurf +rSideSurf;
    addOuterUnionSurf(Out);
    System.addCell(MonteCarlo::Qhull(cellIndex++, fastRefMat, 0.0, Out));
    
    */

    Out = ModelSupport::getComposite(SMap, SI, " -200 -15 ") + tSurf;
    addOuterUnionSurf(Out);
    System.addCell(MonteCarlo::Qhull(cellIndex++, WaterCylinderMat, 0.0, Out));

    return; 
  }

  void
  RecModerator::createLinks()
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

  void RecModerator::createAll(Simulation& System, const attachSystem::FixedComp& FC)
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

    ELog::RegMethod RegA("RecModerator","createAll(Simulation &, const attachSystem::FixedComp&)");
    ELog::EM << "!!! This method should not be called !!!" << ELog::endErr;

    return;
  }



  void RecModerator::createAll(Simulation& System, const attachSystem::FixedComp& FC, const attachSystem::FixedComp& Target, const long int tIndex, const attachSystem::FixedComp& ShutterBay, const long int sIndex, const attachSystem::FixedComp& Reflector, const long int rSide, const long int rTop)
  {
    /*!
      Extrenal build everything
      \param System :: Simulation
      \param FC :: FixedComponent for origin

      In our case the reflector has to be built relative to an origin and an axes set. 
      If you take a simple fixed object, then the axes is the axes set of this fixed object and the origin is the origin of this fixed object,
      which does not mean that the reflector and the object have the same origin, that's just the way you start and then you add the next bits.
    */

    ELog::RegMethod RegA("RecModerator","createAll");
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
  RecModerator::getSurfacePoint(const size_t layerIndex,
				const size_t sideIndex) const
  /*!
    Given a side and a layer calculate the link point
    \param sideIndex :: Side [0-5]
    \param layerIndex :: layer, 0 is inner moderator [0-6]
    \return Surface point
  */
  {
    ELog::RegMethod RegA("RecModerator","getSurfacePoint");

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
  RecModerator::getLayerString(const size_t layerIndex,
			       const size_t sideIndex) const
  /*!
    Given a side and a layer calculate the link surf
    \param sideIndex :: Side [0-5]
    \param layerIndex :: layer, 0 is inner moderator [0-4]
    \return Surface string
  */
  {
    ELog::RegMethod RegA("RecModerator","getLayerString");

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
  RecModerator::getLayerSurf(const size_t layerIndex,
			     const size_t sideIndex) const
  /*!
    Given a side and a layer calculate the link surf
    \param sideIndex :: Side [0-5]
    \param layerIndex :: layer, 0 is inner moderator [0-4]
    \return Surface string
  */
  {
    ELog::RegMethod RegA("RecModerator","getLayerSurf");

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
  RecModerator::getCommonSurf(const size_t sideIndex) const
  /*!
    Given a side calculate the boundary surface
    \param sideIndex :: Side [0-5]
    \return Common dividing surface [outward pointing]
  */
  {
    ELog::RegMethod RegA("RecModerator","getCommonSurf");

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

  RecModerator*
  RecModerator::clone() const
  /*!
    Clone copy constructor
    \return copy of this
  */
  {
    return new RecModerator(*this);
  }


}  // NAMESPACE instrumentSystem
