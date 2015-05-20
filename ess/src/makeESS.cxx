#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <utility>
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
#include <boost/format.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include <boost/range.hpp>

#include "Triple.h"
#include "Exception.h"
#include "FileReport.h"
#include "NameStack.h"
#include "RegMethod.h"
#include "GTKreport.h"
#include "OutputLog.h"
#include "BaseVisit.h"
#include "BaseModVisit.h"
#include "MatrixBase.h"
#include "Matrix.h"
#include "Vec3D.h"
#include "inputParam.h"
#include "Surface.h"
#include "surfIndex.h"
#include "surfRegister.h"
#include "objectRegister.h"
#include "Rules.h"
#include "Code.h"
#include "varList.h"
#include "FuncDataBase.h"
#include "HeadRule.h"
#include "Object.h"
#include "Simulation.h"
#include "LinkUnit.h"
#include "FixedComp.h"
#include "ContainedComp.h"
#include "ContainedGroup.h"
#include "LayerComp.h"
#include "World.h"
#include "FlightLine.h"
#include "AttachSupport.h"
#include "pipeUnit.h"
#include "PipeLine.h"

#include "ProtonVoid.h"
#include "ProtonBeamWindow.h"
#include "WheelBase.h"
#include "Wheel.h"
#include "SegWheel.h"
#include "BeRef.h"
#include "CylinderCell.h"
#include "Grooving.h"
#include "FlowGuide.h"
#include "F5Calc.h"
#include "F5Collimator.h"
#include "CellMap.h"
#include "ModBase.h"
#include "CylModerator.h"
#include "ThermalModerator.h"
#include "RecModerator.h"
#include "ButterflyModerator.h"
#include "BlockAddition.h"
#include "CylPreMod.h"
#include "SupplyPipe.h"
#include "BulkModule.h"
#include "ShutterBay.h"
#include "GuideBay.h"
#include "WaterPipe.h"

#include "stringCombine.h" // for StrFunc

#include "pairRange.h"
#include "NRange.h"
#include "Tally.h"
#include "meshTally.h"

//#include "ConicModerator.h"
#include "essDBMaterial.h"

#include "makeESS.h"

namespace essSystem
{

  makeESS::makeESS() :
    TopReflector(new BeRef("TopBeRef")),
    LowReflector(new BeRef("LowBeRef")),
    OnionModPipe(new FlowGuide("OnionMod")),
    OnionPrePipe(new FlowGuide("OnionPre")),
    PBeam(new essSystem::ProtonVoid("ProtonBeam")),
    PBW(new essSystem::ProtonBeamWindow("ProtonBeamWindow")),
    LowAFL(new essSystem::FlightLine("LowAFlight")),
    LowBFL(new essSystem::FlightLine("LowBFlight")),
    LowPre(new CylPreMod("LowPre")),
    LowSupplyPipe(new constructSystem::SupplyPipe("LSupply")),
    LowReturnPipe(new constructSystem::SupplyPipe("LReturn")),
    TopSupplyPipe(new constructSystem::SupplyPipe("TSupply")),
    TopReturnPipe(new constructSystem::SupplyPipe("TReturn")),

    TopMod(new essSystem::CylModerator("TopMod")), // (new constructSytem::CylMod("TopMod")),
    TopAFL(new essSystem::FlightLine("TopAFlight")),
    TopBFL(new essSystem::FlightLine("TopBFlight")),
    TopPre(new CylPreMod("TopPre")),

    Bulk(new BulkModule("Bulk")),
    //  BulkLowAFL(new essSystem::FlightLine("BulkLAFlight")),
    ShutterBayObj(new ShutterBay("ShutterBay"))
    //    F5(new F5Collimator("F5")),
    //    F15(new F5Collimator("F15"))
    /*!
      Constructor
    */
  {
    ModelSupport::objectRegister& OR=
      ModelSupport::objectRegister::Instance();

    OR.addObject(TopReflector);
    OR.addObject(LowReflector);
    OR.addObject(OnionModPipe);
    OR.addObject(OnionPrePipe);
    OR.addObject(PBeam);
    OR.addObject(PBW);
    OR.addObject(LowAFL);
    OR.addObject(LowBFL);
    OR.addObject(LowPre);
    OR.addObject(TopMod);
    OR.addObject(TopAFL);
    OR.addObject(TopBFL);
    OR.addObject(TopPre);

    OR.addObject(Bulk);
    //  OR.addObject(BulkLowAFL);

    OR.addObject(ShutterBayObj);
  }


  makeESS::~makeESS()
  /*!
    Destructor
  */
  {}

  void makeESS::buildF5Collimator(Simulation& System, std::shared_ptr<F5Collimator> f)
  {
    ELog::RegMethod RegA("makeESS", "buildF5Collimator");
    ModelSupport::objectRegister& OR = ModelSupport::objectRegister::Instance();

    //    std::shared_ptr<F5Collimator> F5Col(new F5Collimator("F5Col"));
    OR.addObject(f);
    f->addInsertCell(74123); // !!! 74123=voidCell tmp - fix it !!!
    f->createAll(System, World::masterOrigin());
    return;

    // find objects by group name
    std::set<int> Oset;
    for (size_t i=0; i<f->NConnect(); i++) {
      const Geometry::Vec3D TestPt = f->getLinkPt(i)*0.99;
      MonteCarlo::Object* OPtr=System.findCell(TestPt, 0);
      if (OPtr) {
	Oset.insert(OPtr->getName());
	ELog::EM << "cell " << OPtr->getName() << " at " << TestPt.X() << " " << TestPt.Y() << " " << TestPt.Z() << ELog::endDiag;
      } else {
	ELog::EM << "no cell at " << TestPt.X() << " " << TestPt.Y() << " " << TestPt.Z() << ELog::endDiag;
      }
    }
    //    return;

    std::set<std::string> ONameSet;
    std::set<int>::const_iterator ac;
    for(ac=Oset.begin();ac!=Oset.end();ac++) {
      ONameSet.insert(OR.inRange(*ac)); // f
      std::cout << "a " << OR.inRange(*ac) << std::endl;
    }

    std::set<std::string>::const_iterator nc;
    for(nc=ONameSet.begin();nc!=ONameSet.end();nc++) {
	const attachSystem::FixedComp* foundObj = OR.getObject<attachSystem::FixedComp>(*nc);
	attachSystem::addToInsertLineCtrl(System, *foundObj, *f);
	std::cout << "c " << *nc << std::endl;
      }
    return;
  }
  
  void
  makeESS::lowFlightLines(Simulation& System)
  /*!
    Build the flight lines of the reflector
    \param System :: Simulation to add to
  */
  {
    ELog::RegMethod RegA("makeESS","lowFlightLines");

    std::string Out;
    Out=LowReflector->getLinkComplement(1)+LowPre->getBoxCut('A');
    LowAFL->addBoundarySurf("liner",Out);  

    Out=LowReflector->getLinkComplement(3)+LowPre->getBoxCut('A');
    LowAFL->addBoundarySurf("linerAl",Out);  

    Out=LowReflector->getLinkComplement(0)+LowPre->getBoxCut('A');
    LowAFL->addBoundarySurf("inner",Out);  
    LowAFL->addBoundarySurf("outer",Out);  
    LowAFL->addOuterSurf("outer",LowPre->getBoxCut('A'));  
    LowAFL->createAll(System,1,*LowPre);
    attachSystem::addToInsertSurfCtrl(System,*LowAFL,*LowPre->getBox('A'));
    attachSystem::addToInsertSurfCtrl(System,*LowReflector, LowAFL->getKey("outer"));

    Out=LowReflector->getLinkComplement(1)+LowPre->getBoxCut('B');
    LowBFL->addBoundarySurf("liner",Out);  

    Out=LowReflector->getLinkComplement(3)+LowPre->getBoxCut('B');
    LowBFL->addBoundarySurf("linerAl",Out);  

    Out=LowReflector->getLinkComplement(0)+LowPre->getBoxCut('B');
    LowBFL->addBoundarySurf("inner",Out);  
    LowBFL->addBoundarySurf("outer",Out);  
    LowBFL->addOuterSurf("outer",LowPre->getBoxCut('B'));  
    LowBFL->createAll(System,0,0,*LowPre);

    attachSystem::addToInsertSurfCtrl(System,*LowBFL,*LowPre->getBox('B'));
    attachSystem::addToInsertSurfCtrl(System,*LowReflector,
				      LowBFL->getKey("outer"));
    //    attachSystem::addToInsertSurfCtrl(System,*LowBFL, LowAFL->getKey("outer"));
    return;
  }

  void
  makeESS::topFlightLines(Simulation& System, const std::string topModType)
  /*!
    Build the flight lines of the reflector
    \param System :: Simulation to add to
  */
  {
    ELog::RegMethod RegA("makeESS","topFlightLines");
    std::string Out;
    Out=TopReflector->getLinkComplement(1)+TopPre->getBoxCut('A');
    TopAFL->addBoundarySurf("liner",Out);  

    Out=TopReflector->getLinkComplement(3)+TopPre->getBoxCut('A');
    TopAFL->addBoundarySurf("linerAl",Out);  

    Out=TopReflector->getLinkComplement(0)+TopPre->getBoxCut('A');
    TopAFL->addBoundarySurf("inner",Out);  
    TopAFL->addBoundarySurf("outer",Out);  
    TopAFL->addOuterSurf("outer",TopPre->getBoxCut('A'));  
    if (topModType=="M3A") {
      TopAFL->createAll(System,0,0,*TopMod);
      attachSystem::addToInsertSurfCtrl(System,*TopPre, TopAFL->getKey("outer"));  // exclude from TopPre
    } else {
      TopAFL->createAll(System,0,1,*TopPre);
    }
    attachSystem::addToInsertSurfCtrl(System,*TopAFL,*TopPre->getBox('A'));
    attachSystem::addToInsertSurfCtrl(System,*TopReflector, TopAFL->getKey("outer"));

    Out=TopReflector->getLinkComplement(1)+TopPre->getBoxCut('B');
    TopBFL->addBoundarySurf("liner",Out);  

    Out=TopReflector->getLinkComplement(3)+TopPre->getBoxCut('B');
    TopBFL->addBoundarySurf("linerAl",Out);  

    Out=TopReflector->getLinkComplement(0)+TopPre->getBoxCut('B');
    TopBFL->addBoundarySurf("inner",Out);  
    TopBFL->addBoundarySurf("outer",Out);  
    TopBFL->addOuterSurf("outer",TopPre->getBoxCut('B'));  
    if (topModType=="M3A") {
      TopBFL->createAll(System,0,1,*TopMod);
      attachSystem::addToInsertSurfCtrl(System,*TopPre, TopBFL->getKey("outer"));  // exclude from TopPre
    } else {
      TopBFL->createAll(System,0,0,*TopPre);
    }
    attachSystem::addToInsertSurfCtrl(System,*TopBFL,*TopPre->getBox('B'));
    attachSystem::addToInsertSurfCtrl(System,*TopReflector,
				      TopBFL->getKey("outer"));

    // Out=TopReflector->getLinkComplement(0)+TopPre->getBoxCut('B');
    // TopBFL->addBoundarySurf("inner",Out);  
    // TopBFL->addBoundarySurf("outer",Out);  
    // TopBFL->createAll(System,0,0,*TopPre);
    // TopBFL->addOuterSurf("outer",TopPre->getBoxCut('B'));  
    // attachSystem::addToInsertSurfCtrl(System,*TopBFL,*TopPre->getBox('B'));
    // attachSystem::addToInsertSurfCtrl(System,*TopReflector,
    // 				    TopBFL->getKey("outer"));


    return;
  }


void 
makeESS::makeTarget(Simulation& System, const mainSystem::inputParam& IParam)
  /*!
    Build the different ESS targets
    \param System :: Simulation
    \param IParam :: Parameters
  */
  {
    ELog::RegMethod RegA("makeESS","makeTarget");

    const int voidCell(74123);  
    const std::string targetType = IParam.getValue<std::string>("targetType");
    if (targetType=="help") {
      ELog::EM<<"Target Type : "<<ELog::endBasic;
      ELog::EM<<"Wheel       : Simple wheel form"<<ELog::endBasic;
      ELog::EM<<"SegWheel    : Segmented wheel"<<ELog::endBasic;
      throw ColErr::ExitAbort("help exit");
    }
    if (targetType=="Wheel")
      Target=std::shared_ptr<WheelBase>(new Wheel("Wheel"));
    else if (targetType=="SegWheel")
      Target=std::shared_ptr<WheelBase>(new SegWheel("SegWheel"));
    else
      throw ColErr::InContainerError<std::string>(targetType,"Unknown target type");

    Target->addInsertCell("Shaft",voidCell); // we put shaft in void because its part is in void
    Target->addInsertCell("Wheel",voidCell); // because back side of the wheel in void
    Target->createAll(System,World::masterOrigin());

    return;
  }


  void 
  makeESS::createGuides(Simulation& System, const std::string lowModType)
  /*!
    Create all the guidebays and guides
    \param System :: Simulation system to use
  */
  {
    ELog::RegMethod RegA("makeESS","createGuides");

    for(size_t i=0;i<4;i++)
      {
	std::shared_ptr<GuideBay> GB(new GuideBay("GuideBay",i+1));
	GB->addInsertCell("Inner",ShutterBayObj->getMainCell());
	GB->addInsertCell("Outer",ShutterBayObj->getMainCell());
	GB->setCylBoundary(Bulk->getLinkSurf(2),
			   ShutterBayObj->getLinkSurf(2));
	if(i<2) {
	  //	  std::cout << "*** crush here with only one moderator  " << std::endl;
	  //	  if (lowModType != "None")	    GB->createAll(System,*LowMod);  
	} else {
	  GB->createAll(System,*TopMod);
	}
	GBArray.push_back(GB);
      }
    return;
  }

  void
  makeESS::buildLowMod(Simulation& System)
  /*!
    Build the lower moderators
    \param System :: Simulation to build
  */
  {
    ELog::RegMethod RegA("makeESS","buildLowMod");

    ModelSupport::objectRegister& OR=
      ModelSupport::objectRegister::Instance();
    LowMod=std::shared_ptr<essSystem::CylModerator>(new essSystem::CylModerator("LowMod"));
    OR.addObject(LowMod);

    LowReflector->addToInsertChain(*LowMod);
    //    LowMod->addInsertCell(LowReflector->getMainCell());
    LowMod->createAll(System, World::masterOrigin()); //*LowReflector);

    // ??? why not use LowReflector->addToInsertChain(LowPre->getKey("Main")) ??? !!!
    for (int i=LowReflector->getRefIndex()+1; i<LowReflector->getCellIndex(); i++) {
      LowPre->addInsertCell("Main", i);
      LowPre->addInsertCell("BlockA", i);
      LowPre->addInsertCell("BlockB", i);
    }

    /*    LowPre->addInsertCell("Main", LowReflector->getMainCell());
    LowPre->addInsertCell("Main", LowReflector->getWallCell());
    LowPre->addInsertCell("BlockA",LowReflector->getMainCell());
    LowPre->addInsertCell("BlockB",LowReflector->getMainCell());*/
    LowPre->createAll(System,*LowMod, *Target, 2);
    //    LowPre->addToInsertChain(Target->getKey("Wheel")); // because target wheel might cut through the gap in Be - this avoids geometry errors

    return;
  }

  void 
  makeESS::dumpMaterialMesh(Simulation& SimPtr, const Geometry::Vec3D &startPt, const Geometry::Vec3D &endPt,
			    const size_t nX, const size_t nY, const size_t nZ,
			    const char *fname)
  {
    /*
      Dumps a mesh with materials in ASCII file 'fname'
      \param startPt :: Point with min coordinates
      \param endPt :: Point with max coordinates
      \param nX :: x-coordinate division
      \param nY :: y-coordinate division
      \param nZ :: z-coordinate division
      \param fname :: output file name
    */
    
    //    nPts

    Geometry::Vec3D Origin = startPt; // start corner (x,y,z=min)
    Geometry::Vec3D XYZ = endPt-Origin;
    Geometry::Vec3D aVec;
    MonteCarlo::Object *ObjPtr(0);
    
    Triple<long int> nPts = Triple<long int>(static_cast<long int>(nX), static_cast<long int>(nY), static_cast<long int>(nZ));
    double  stepXYZ[3];
    for (size_t i=0; i<3; i++)
      stepXYZ[i] = XYZ[i]/static_cast<double>(nPts[i]);
  

    std::ofstream fmesh;
    fmesh.open(fname);
    int mat = 0;

    for (size_t i=0; i<nX; i++) {
      aVec[0] = stepXYZ[0]*(static_cast<double>(i)+0.5);
      for (size_t j=0; j<nY; j++) {
	aVec[1] = stepXYZ[1]*(static_cast<double>(j)+0.5);
	for (size_t k=0; k<nZ; k++) {
	  aVec[2] = stepXYZ[2]*(static_cast<double>(k)+0.5);
	  const Geometry::Vec3D Pt = Origin+aVec;
	  ObjPtr = SimPtr.findCell(Pt, ObjPtr);
	  if (ObjPtr)
	    mat = ObjPtr->getMat();
	  else 
	    mat = -1;
	  fmesh << Pt[0] << " " << Pt[1] << " " << Pt[2] << "\t" << mat << std::endl;
	  
	}
      }
    }
    fmesh.close();
    
  }


  void 
  makeESS::build(Simulation* SimPtr,
		 const mainSystem::inputParam& IParam)
  /*!
    Carry out the full build
    \param SimPtr :: Simulation system
    \param IParam :: Input parameters
  */
  {
    ELog::RegMethod RControl("makeESS","build");    // For output stream


    const std::string lowModType=IParam.getValue<std::string>("lowMod"); // -lowMod
    const std::string topModType=IParam.getValue<std::string>("topMod");

    const std::string topModCooling=IParam.getValue<std::string>("topModCooling");
    const std::string lowModCooling=IParam.getValue<std::string>("lowModCooling");
    const std::string topPreCooling=IParam.getValue<std::string>("topPreCooling");
    const std::string isLowWaterDisc=IParam.getValue<std::string>("lowWaterDisc");
    const std::string isTopWaterDisc=IParam.getValue<std::string>("topWaterDisc");
    const std::string isThermalCylMod=IParam.getValue<std::string>("thermalCylMod");
    const int nF5 = IParam.getValue<int>("nF5");

    ModelSupport::objectRegister& OR = ModelSupport::objectRegister::Instance();

    bool isf5col = true; // set to true to enable F5 collimator

    int voidCell(74123); // MegaVoid - magic number which should not be changed
    // Add extra materials to the DBdatabase
    ModelSupport::addESSMaterial();

    makeTarget(*SimPtr,IParam);
  
    TopReflector->addInsertCell(voidCell);
    TopReflector->createAll(*SimPtr,World::masterOrigin(), *Target, 3);
    TopReflector->addToInsertChain(Target->getKey("Wheel")); // because target wheel cuts through reflector

    Bulk->createAll(*SimPtr, World::masterOrigin(), *TopReflector, *LowReflector);
    Bulk->addToInsertChain(*TopReflector);
    Bulk->addToInsertChain(*LowReflector);

    attachSystem::addToInsertSurfCtrl(*SimPtr,*Bulk,Target->getKey("Wheel")); // [2:865] test every individual cell that makes up the bulk to see which one is intersect wheel and the do the intersection if it is required. it's slow, so should not be used when no needed
    attachSystem::addToInsertForced(*SimPtr,*Bulk,Target->getKey("Shaft"));

    LowReflector->addInsertCell(voidCell);
    LowReflector->createAll(*SimPtr, World::masterOrigin(), *Target, 2);
    LowReflector->addToInsertChain(Target->getKey("Wheel")); // because target wheel cuts through reflector


    if (lowModType=="Cone") {
      //buildConicMod(*SimPtr);
      //Bulk->addFlightUnit(*SimPtr,*LowAFL, *LowReflector);
      //Bulk->addFlightUnit(*SimPtr,*LowBFL, *LowReflector);
    } else if ((lowModType=="Tube") or (lowModType=="GroovedTube")) {  // Tube moderator is a thermal moderator with complicated water structure
      TubeMod = std::shared_ptr<RecModerator>(new RecModerator("TubeMod"));
      OR.addObject(TubeMod);
      LowReflector->addToInsertChain(*TubeMod);
      Bulk->addToInsertChain(*TubeMod);
      TubeMod->createAll(*SimPtr, *LowReflector, *Target, 2, *Bulk, -3, *LowReflector, -2, -5);

      if (lowModType=="GroovedTube") {
	//	std::cout << "Here " << lowModType << std::endl;
      // grooving
	TubeModGrooving = std::shared_ptr<Grooving>(new Grooving("Groove"));
	OR.addObject(TubeModGrooving);
	TubeMod->addToInsertChain(*TubeModGrooving);
	LowReflector->addToInsertChain(*TubeModGrooving);
	TubeModGrooving->createAll(*SimPtr, *TubeMod, *TubeMod, 0, 1);
      }

    } else if (lowModType == "Butterfly") {
      if ((isLowWaterDisc == "On") && (!LowWaterDisc)) { // requested but not yet built
	LowWaterDisc = std::shared_ptr<BeRef>(new BeRef("LowWaterDisc"));
	OR.addObject(LowWaterDisc);
	LowReflector->addToInsertChain(*LowWaterDisc);
	LowWaterDisc->createAll(*SimPtr, *LowReflector, *Target, 2);
      }

      LowButterfly = std::shared_ptr<ButterflyModerator>(new ButterflyModerator("LowButterfly"));
      OR.addObject(LowButterfly);
      LowReflector->addToInsertChain(*LowButterfly); // !!! kbat
      Bulk->addToInsertChain(*LowButterfly);
      if (LowWaterDisc) LowWaterDisc->addToInsertChain(*LowButterfly);
      LowButterfly->setReflectorSurfaces(*LowReflector, -2, -4); // -2 -4
      LowButterfly->createAll(*SimPtr, *LowReflector, *Bulk, -3); // -3

      if ((LowButterfly->IsTopPreCoolingChannels()) && (LowButterfly->getTopPreType()==2)) { // cooling channels implemented for TopPreType==2 only
	LowModOnion = std::shared_ptr<FlowGuide>(new FlowGuide("LowModOnion"));
	OR.addObject(LowModOnion);
	LowModOnion->setUpperSurface(*LowButterfly, -9);
	LowModOnion->setBottomSurface(*LowButterfly, 5);
	LowButterfly->addToInsertChain(*LowModOnion);
	LowModOnion->createAll(*SimPtr, *LowButterfly);
      }

      if (LowButterfly->getPipeType()==1) {
	LowSupplyPipe->createAll(*SimPtr,*LowButterfly,0,2,1, *LowButterfly, 1);
	LowReturnPipe->createAll(*SimPtr,*LowButterfly,0,1,0, *LowButterfly, 0);
      }

    } else if (lowModType == "Thermal") { // thermal box
      ThermalRecMod = std::shared_ptr<ThermalModerator>(new ThermalModerator("ThermalRecMod"));
      OR.addObject(ThermalRecMod);
      LowReflector->addToInsertChain(*ThermalRecMod);
      Bulk->addToInsertChain(*ThermalRecMod);
      //      ThermalRecMod->createAll(*SimPtr, *LowReflector, *Target, 2, *Bulk, -3, *LowReflector, -2, -5);
      ThermalRecMod->createAll(*SimPtr, *LowReflector, *LowReflector, 5, *Bulk, -3, *LowReflector, -2, -5);
    } else if (lowModType=="None") {
    }   else { // cylindrical = "Basic"
      buildLowMod(*SimPtr);
      if (isLowWaterDisc == "On") {
	LowWaterDisc = std::shared_ptr<BeRef>(new BeRef("LowWaterDisc"));
	OR.addObject(LowWaterDisc);
	LowReflector->addToInsertChain(*LowWaterDisc);
	LowPre->addToInsertChain(*LowWaterDisc);
	LowMod->addToInsertChain(*LowWaterDisc);
	LowWaterDisc->addInsertCell(voidCell); // for Low supply/return pipes
	LowWaterDisc->createAll(*SimPtr, *LowReflector, *Target, 2);
      }
      if (isThermalCylMod == "On") {
	ThermalCylMod = std::shared_ptr<CylinderCell>(new CylinderCell("ThermalCylMod"));
	OR.addObject(ThermalCylMod);
	LowReflector->addToInsertChain(*ThermalCylMod);
	LowPre->addToInsertChain(*ThermalCylMod);
	LowMod->addToInsertChain(*ThermalCylMod);
	ThermalCylMod->createAll(*SimPtr, *LowReflector);
      }
      lowFlightLines(*SimPtr);
      Bulk->addFlightUnit(*SimPtr,*LowAFL, *LowReflector);
      Bulk->addFlightUnit(*SimPtr,*LowBFL, *LowReflector); 
    }


    if (topModType == "Butterfly") {
      if ((isTopWaterDisc == "On") && (!TopWaterDisc)) { // requested but not yet built
	TopWaterDisc = std::shared_ptr<BeRef>(new BeRef("TopWaterDisc"));
	OR.addObject(TopWaterDisc);
	TopReflector->addToInsertChain(*TopWaterDisc);
	TopWaterDisc->createAll(*SimPtr, *TopReflector, *Target, 3);
      }

      TopButterfly = std::shared_ptr<ButterflyModerator>(new ButterflyModerator("TopButterfly"));
      OR.addObject(TopButterfly);
      TopReflector->addToInsertChain(*TopButterfly); // !!! kbat
      Bulk->addToInsertChain(*TopButterfly);
      if (TopWaterDisc) TopWaterDisc->addToInsertChain(*TopButterfly);

      TopButterfly->setReflectorSurfaces(*TopReflector, -2, -4);
      TopButterfly->createAll(*SimPtr, *TopReflector, *Bulk, -3);

      if (TopButterfly->getPipeType()==1) { // horisontal
	TopSupplyPipe->addInsertCell(4, TopWaterDisc->getCell("TopWaterDiscRing",1)); // 4=pipe layer number, 1=ring number
	TopSupplyPipe->addInsertCell(4, TopButterfly->getCell("TopButterflyWall", 1));
	TopSupplyPipe->createAll(*SimPtr,*TopButterfly,0,2,1, *TopButterfly, 1);
	
	TopReturnPipe->addInsertCell(4, TopWaterDisc->getCell("TopWaterDiscRing",1)); // 4=pipe layer number, 1=ring number
	TopReturnPipe->addInsertCell(4, TopButterfly->getCell("TopButterflyWall", 1));
	TopReturnPipe->createAll(*SimPtr,*TopButterfly,0,1,0, *TopButterfly, 0);
      } else if (TopButterfly->getPipeType()==2) { // vertical
	TopSupplyPipe->createAll(*SimPtr,*TopButterfly,0,6,5, *TopButterfly, 5);
	TopReturnPipe->createAll(*SimPtr,*TopButterfly,0,6,5, *TopButterfly, 5);
      }

    } else {
      TopReflector->addToInsertChain(*TopMod);
      TopMod->createAll(*SimPtr, World::masterOrigin()); //*TopReflector);//, *Target);//->getKey("Wheel"));
      for (int i=TopReflector->getRefIndex()+1; i<TopReflector->getCellIndex(); i++) {
	TopPre->addInsertCell("Main", i);
	TopPre->addInsertCell("BlockA", i);
	TopPre->addInsertCell("BlockB", i);
      }

      TopPre->createAll(*SimPtr,*TopMod, *Target, 3);
    //    TopPre->addToInsertChain(Target->getKey("Wheel")); // because target wheel might cut through the gap in Be
      topFlightLines(*SimPtr, topModType);

      Bulk->addFlightUnit(*SimPtr,*TopAFL, *TopReflector);
      Bulk->addFlightUnit(*SimPtr,*TopBFL, *TopReflector);
    }


    // Full surround object
    ShutterBayObj->addInsertCell(voidCell);
    ShutterBayObj->createAll(*SimPtr,*Bulk,*Bulk);
    attachSystem::addToInsertForced(*SimPtr,*ShutterBayObj,
				    Target->getKey("Wheel"));
    attachSystem::addToInsertForced(*SimPtr,*ShutterBayObj,
				    Target->getKey("Shaft"));

    TopReflector->addToInsertChain(*PBeam); // because proton beam cuts through reflector
    LowReflector->addToInsertChain(*PBeam); // because proton beam cuts through reflector
    Bulk->addToInsertChain(*PBeam); // because proton beam cuts through bulk
    ShutterBayObj->addToInsertChain(*PBeam);
    PBeam->createAll(*SimPtr,*Target,1,*ShutterBayObj,-3); // prolong pbeam up to shutter bay outer surface
  
    PBeam->addToInsertChain(*PBW);
    ShutterBayObj->addToInsertChain(*PBW);
    PBW->createAll(*SimPtr, *PBeam, -1);

  
 
    if (topModCooling == "Onion") {
      TopMod->addToInsertChain(*OnionModPipe);
      OnionModPipe->setUpperSurface(*TopMod, -8);
      OnionModPipe->setBottomSurface(*TopMod, -7);
      OnionModPipe->createAll(*SimPtr, *TopMod);
      if (topModType == "M3A") {
	//	std::cout << "hereA" << std::endl;
	TopSupplyPipe->createAll(*SimPtr,*TopMod,0,4,3, *TopMod,3);// works: 043,3
	// orig: TopReturnPipe->createAll(*SimPtr,*TopMod,0,3,2, *TopMod,2); // works: 032,2
	// play here:
	TopReturnPipe->createAll(*SimPtr,*TopMod,0,3,2, *TopMod, 2);
      } else {
	TopSupplyPipe->createAll(*SimPtr,*TopMod,0,4,3, *TopPre,5); // works: 0,4,3 5
	TopReturnPipe->createAll(*SimPtr,*TopMod,0,3,2, *TopPre,4); // works: 0,3,2 4
      }

    } else if (topModCooling=="Standard") {
      TopSupplyPipe->createAll(*SimPtr,*TopMod,0,6,4, *TopPre,4); // 0644 is ok for vertical
      //TopReturnPipe->createAll(*SimPtr,*TopMod,0,3,2, *TopPre,4);
    }

    if (lowModCooling=="Standard") {
      LowSupplyPipe->createAll(*SimPtr,*LowMod,0,6,4, *LowPre,2);
    }

    // PreModerator onion cooling
    if (topPreCooling=="Onion") {
      TopPre->addToInsertChain(*OnionPrePipe);
      OnionPrePipe->setUpperSurface(*TopPre, -8);
      OnionPrePipe->setBottomSurface(*TopPre, -7);
      TopReflector->addToInsertChain(*OnionPrePipe); // just in case PrePipe goes into reflector
      OnionPrePipe->createAll(*SimPtr, *TopPre);
      
      // since both OnionPrePipe and TopMod have already been created, we can't
      // call OnionPrePipe::addToInsertChain anymore => use addToInsertSurfCtrl:
      attachSystem::addToInsertSurfCtrl(*SimPtr,  *OnionPrePipe, *TopMod); // onion pre pipes do not belong to TopMod


      // WATER PIPES (one object for supply/return)

      TopPreAPipe = std::shared_ptr<WaterPipe>(new WaterPipe("TPreA"));
      TopPreBPipe = std::shared_ptr<WaterPipe>(new WaterPipe("TPreB"));

      OR.addObject(TopPreAPipe);
      OR.addObject(TopPreBPipe);

      TopPreAPipe->createAll(*SimPtr, *TopAFL, 0, *TopAFL, 0);
      TopPreBPipe->createAll(*SimPtr, *TopBFL, 0, *TopBFL, 0);
      
    }

    //createGuides(*SimPtr, lowModType);
    // F5 collimator
    if (isf5col) { 
      for (int i=0; i<nF5; i++) {
	std::shared_ptr<F5Collimator> F5(new F5Collimator(StrFunc::makeString("F", i*10+5).c_str()));
	buildF5Collimator(*SimPtr, F5);
	attachSystem::addToInsertSurfCtrl(*SimPtr, *ShutterBayObj, *F5);
	F5array.push_back(F5);
      }
    }


    // write out the material mesh
    if (IParam.flag("matmesh"))
      {
	Geometry::Vec3D ptStart(-40, -40, -20);
	Geometry::Vec3D ptEnd(40, 40, 20);
        dumpMaterialMesh(*SimPtr, ptStart, ptEnd, 200, 200, 100, "mesh.dat");
      }
      
  }

}   // NAMESPACE essSystem

