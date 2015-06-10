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
#include "ButterflyModerator.h"
#include "Plane.h"

#include <boost/algorithm/string/replace.hpp>

namespace essSystem
{

  ButterflyModerator::ButterflyModerator(const std::string& Key) :
    constructSystem::ModBase(Key, 9)
    /*!
      Constructor
      \param Key :: Name of construction key
    */
  {}

  ButterflyModerator::ButterflyModerator(const ButterflyModerator& A) : 
    ModBase(A)
    /*!
      Copy constructor
      \param A :: ButterflyModerator to copy
    */
  {}



  ButterflyModerator::~ButterflyModerator()
  /*!
    Destructor
  */
  {}

  void
  ButterflyModerator::populate(const FuncDataBase& Control)
  /*!
    Populate all the variables
    \param Control :: Variable table to use
  */
  {
    ELog::RegMethod RegA("ButterflyModerator","populate");

    // Master values
    xStep=Control.EvalVar<double>(keyName+"XStep");
    yStep=Control.EvalVar<double>(keyName+"YStep");
    zStep=Control.EvalVar<double>(keyName+"ZStep");
    xyAngle=Control.EvalDefVar<double>(keyName+"XYangle", 0.0);
    zAngle=Control.EvalDefVar<double>(keyName+"Zangle", 0.0);
    nLayers=Control.EvalVar<size_t>(keyName+"NLayers");   

    double thick, W, D, H, T;
    int M;
    for (size_t i=0; i<nLayers; i++) {
      if (i) {
	// *2 because container from both sides
	thick = 2*Control.EvalDefVar<double>(StrFunc::makeString(keyName+"Thick", i), 0.2);
	//	D += 2*Control.EvalDefVar<double>(StrFunc::makeString(keyName+"YGap", i), 0.2);
	//	H += 2*Control.EvalDefVar<double>(StrFunc::makeString(keyName+"ZGap", i), 0.2);
	M = ModelSupport::EvalDefMat<int>(Control, StrFunc::makeString(keyName+"Mat", i), 0);
	T = (!M) ? 0.0 : Control.EvalDefVar<double>(StrFunc::makeString(keyName+"Temp", i), 0.0);
	
	W += thick;
	D += thick;
	H += thick;

	Width.push_back(W);
	Length.push_back(D);
	Height.push_back(H);
      } else {
	W = Control.EvalDefVar<double>(keyName+"Width", 10);
	D = Control.EvalDefVar<double>(keyName+"Length", 10);
	H = Control.EvalDefVar<double>(keyName+"Height", 3);
	M = ModelSupport::EvalDefMat<int>(Control, keyName+"ModMat", 0);
	T = Control.EvalDefVar<double>(keyName+"Temp", 0);

	Width.push_back(W);
	Length.push_back(D);
	Height.push_back(H);
      }
      mat.push_back(M);
      temp.push_back(T);
    }


    WingLength = Control.EvalVar<double>(keyName+"WingLength");
    WingCurvatureSideRadius = Control.EvalVar<double>(keyName+"WingCurvatureSideRadius");
    WingCurvatureCentralRadius = Control.EvalVar<double>(keyName+"WingCurvatureCentralRadius");

    PreMat = ModelSupport::EvalMat<int>(Control, keyName+"PreMat");
    PreTemp = Control.EvalVar<double>(keyName+"PreTemp"); 
    PreLength = Control.EvalVar<double>(keyName+"PreLength"); 
    PreThick1 = Control.EvalVar<double>(keyName+"PreThick1");
    PreThick2 = Control.EvalVar<double>(keyName+"PreThick2");
    PreBaseThick = Control.EvalVar<double>(keyName+"PreBaseThick");
    PreHeight = Control.EvalVar<double>(keyName+"PreHeight");
    PreWallThick = Control.EvalVar<double>(keyName+"PreWallThick");
    PreWallMat = Control.EvalVar<int>(keyName+"PreWallMat");
    PreWallTemp = Control.EvalVar<double>(keyName+"PreWallTemp");

    PreTopHeight = Control.EvalVar<double>(keyName+"PreTopHeight");
    PreBottomHeight = Control.EvalVar<double>(keyName+"PreBottomHeight");
    

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
  
    FlightXOffset = Control.EvalVar<double>(keyName+"FlightXOffset");
    FlightLineType = Control.EvalDefVar<int>(keyName+"FlightLineType", 0);
    if (FlightLineType>1) ELog::EM << "FlightLineType must be <= 1: " << FlightLineType << ELog::endErr;
    FlightLineWrapTopPreType = Control.EvalDefVar<int>(keyName+"FlightLineWrapTopPreType", 0);
    if (FlightLineWrapTopPreType>1) ELog::EM << "FlightLineWrapTopPreType must be <= 1: " << FlightLineWrapTopPreType << ELog::endErr;

    TopPreType = Control.EvalDefVar<int>(keyName+"TopPreType", 0.0);
    if (TopPreType>3)
      ELog::EM << "TopPreType must be <= 3: " << TopPreType << ELog::endErr;
    TopPreRadius = Control.EvalDefVar<double>(keyName+"TopPreRadius", 39.0/2);
    TopPreHeight = Control.EvalDefVar<double>(keyName+"TopPreHeight", 1);
    TopPreWidth = Control.EvalDefVar<double>(keyName+"TopPreWidth", 40);

    TopPreCoolingChannels = Control.EvalDefVar<int>(keyName+"TopPreCoolingChannels", 1);
    TopPreTopVoidHeight = Control.EvalDefVar<double>(keyName+"TopPreTopVoidHeight", 0.0);
    TopPreTopAlHeight = Control.EvalDefVar<double>(keyName+"TopPreTopAlHeight", 0.0);

    MaterialInBeRef = ModelSupport::EvalDefMat<int>(Control, keyName + "MaterialInBeRef", 1011);

    PipeType = Control.EvalVar<int>(keyName+"PipeType");

    return;
  }

  void
  ButterflyModerator::createUnitVector(const attachSystem::FixedComp& FC)
  /*!
    Create the unit vectors
    \param FC :: Fixed Component
  */
  {
    ELog::RegMethod RegA("ButterflyModerator","createUnitVector");
    attachSystem::FixedComp::createUnitVector(FC);
    applyShift(xStep,yStep,zStep);
    applyAngleRotate(xyAngle,zAngle);
  
    return;
  }

  void ButterflyModerator::setReflectorSurfaces(const attachSystem::FixedComp& Reflector, const long int rSideAl, const long int rSideBe)
  {
    ReflectorSideAl = (rSideAl<0) ? Reflector.getLinkComplement(static_cast<size_t>(-(rSideAl+1))) : Reflector.getLinkString(static_cast<size_t>(rSideAl));
    ReflectorSideBe = (rSideBe<0) ? Reflector.getLinkComplement(static_cast<size_t>(-(rSideBe+1))) : Reflector.getLinkString(static_cast<size_t>(rSideBe));
  }




  /*  int sign(double f)
  {
    return (f<0) ? -1 : 1;
    }*/

  void
  ButterflyModerator::createSurfaces() 
  {
    /*
      // see archetypes/butterfly.svgz

      x+y- 30      20 x+y+
             _      _
            / \    / \
            |  8 7   |
            3   -    4
            |  2 1   |
            \_/   \_/
       x-y- 40     10 x-y+
     */
    ELog::RegMethod RegA("ButterflyModerator","createSurfaces");

    //    std::ofstream cylph;
    //    cylph.open("set_variables.dat", std::ios::out);

    double WingRad = WingCurvatureSideRadius; // wing curvature radius in the given layer
    std::vector<double> vWingRad; // array of radii - used to round corners of top premoderator
    double WingCentralRad = WingCurvatureCentralRadius;
    double dRad = 0.0; // additional shift in WingRad and coordinates
    
    double fullWidth = 0.0;
    double theta = 0.0; // angle of inclination of the inclined planes in the first layer (to make equal thickness of the other layers )
    struct {
      double x;
      double y;
    } A[8], B[4], C[2], D[6], E[2], F[4], G[4], K[4];
    Geometry::Vec3D dirX(X);
    Geometry::Plane *sideWingPlane[4]; // side wing planes
    double virtAngle[4]; // angles of virtual planes with respect to the positive y-axis
    double sideAngle[4]; // angles of side planes 1, 2, 7, 8

    int SI(modIndex);
    for (size_t i=0; i<nLayers; i++) {
      if (i>0) {
	dRad = (Width[i]-Width[i-1])/2;
	WingRad += dRad;
	WingCentralRad += dRad;
      }
      vWingRad.push_back(WingRad);
      double xangle = atan(Width[i]/2/WingLength);
      std::cout << "xangle: " << xangle*180/M_PI << std::endl;

      ModelSupport::buildPlane(SMap, SI+3, Origin-Y*(Length[i]/2.0), Y);
      ModelSupport::buildPlane(SMap, SI+4, Origin+Y*(Length[i]/2.0), Y); // this one at Y<0
      ModelSupport::buildPlane(SMap, SI+5, Origin-Z*(Height[i]/2.0), Z);
      ModelSupport::buildPlane(SMap, SI+6, Origin+Z*(Height[i]/2.0), Z);

      // todo: create these planes before the for... loop
      ModelSupport::buildPlane(SMap, SI+11, Origin, Y); // plane Y=0
      ModelSupport::buildPlane(SMap, SI+12, Origin, Z); // plane Z=0

      // peak of the right wing - where inclined planes intersept
      C[0].x = 0; // center
      C[0].y = Length[i]/2-WingLength;
      if (i>0) {
	theta = 0.0; //xangle;
	C[0].y -= (Width[i]-Width[0])/sin(xangle); // here must be sin (checked)
      }

      // peak of the left wing
      C[1].x = C[0].x;
      C[1].y = -C[0].y; // total symmetry :)

      // ROUNDING THE CORNERS
      // cylinder in the x-y+ corner
      D[0].x = -Width[0]/2.0+WingCurvatureSideRadius;
      D[0].y = Length[0]/2.0-WingCurvatureSideRadius;
      ModelSupport::buildCylinder(SMap, SI+9, Origin+X*(D[0].x)-Y*(D[0].y), Z, WingRad);
      
      // intersect of cylinder 9 and plane 1 !!! todo make it more accurate - find intersept of a tangential !!!
      A[0].x = D[0].x-WingRad*cos(xangle);
      A[0].y = D[0].y-WingRad*sin(xangle);

      // intersect of cylinder 9 and plane 4
      B[0].x = D[0].x+0;
      B[0].y = D[0].y+WingRad;

      D[1].x = Width[i]/2-WingRad;
      D[1].y = D[0].y;

      /*      std::cout << "A0 " << A[0].x << " " << A[0].y << std::endl;
      std::cout << "C0 " << C[0].x << " " << C[0].y << std::endl;
      std::cout << "D0 " << D[0].x << " " << D[0].y << std::endl;*/


      dirX = X;
      sideAngle[0] = atan2(A[0].x-C[0].x, A[0].y-C[0].y)*180/M_PI;
      Geometry::Quaternion::calcQRotDeg(sideAngle[0], Z).rotate(dirX);
      sideWingPlane[0] = ModelSupport::buildPlane(SMap, SI+1, Origin+X*A[0].x - Y*A[0].y, dirX);

      dirX = X;
      virtAngle[0] = atan2(B[0].x-A[0].x, B[0].y-A[0].y)*180/M_PI;
      Geometry::Quaternion::calcQRotDeg(virtAngle[0], Z).rotate(dirX);
      ModelSupport::buildPlane(SMap, SI+10, Origin+X*A[0].x - Y*A[0].y, dirX);
      

      // cylinder in the x+y+ corner
      ModelSupport::buildCylinder(SMap, SI+19, Origin+X*(D[1].x)-Y*(D[1].y), Z, WingRad);

      // intersect of cylinder 19 and plane 7
      A[1].x = D[1].x+WingRad*cos(xangle);
      A[1].y = D[1].y-WingRad*sin(xangle);
      //      std::cout << "A1 " << Origin+X*A[1].x - Y*A[1].y << std::endl; // ok
      
      // intersect of cylinder 19 and plane 4
      B[1].x = D[1].x;
      B[1].y = D[1].y+WingRad;
      //      std::cout << "B1 " << Origin+X*B[1].x - Y*B[1].y << std::endl; // ok

      dirX = X;
      sideAngle[1] = atan2(A[1].x-C[0].x, A[1].y-C[0].y)*180/M_PI;
      Geometry::Quaternion::calcQRotDeg(sideAngle[1], Z).rotate(dirX);
      //      sideWingPlane[1] = ModelSupport::buildPlane(SMap, SI+7, Origin+Geometry::Vec3D(A[1].x, A[1].y, 0), dirX);
      sideWingPlane[1] = ModelSupport::buildPlane(SMap, SI+7, Origin+X*A[1].x - Y*A[1].y, dirX);

      dirX = X;
      virtAngle[1] = atan2(B[1].x-A[1].x, B[1].y-A[1].y)*180/M_PI;
      Geometry::Quaternion::calcQRotDeg(virtAngle[1], Z).rotate(dirX);
      //      ModelSupport::buildPlane(SMap, SI+20, Origin+Geometry::Vec3D(A[1].x, A[1].y, 0), dirX);
      ModelSupport::buildPlane(SMap, SI+20, X*A[1].x -Y*A[1].y, dirX);

      // cylinder in the x+y- corner
      D[2].x = D[1].x;
      D[2].y = -Length[i]/2+WingRad;
      ModelSupport::buildCylinder(SMap, SI+29, Origin+X*(D[2].x)-Y*(D[2].y), Z, WingRad);
      //      std::cout << "c29 " << Origin+X*(D[2].x)-Y*(D[2].y) << std::endl; // OK


      // intersect of cylinder 29 and plane 8
      A[2].x = A[1].x;
      A[2].y = D[2].y+WingRad*sin(xangle);
      //      std::cout << "A2 " << Origin+X*A[2].x - Y*A[2].y << std::endl; // ok

      // intersect of cylinder 19 and plane 3
      B[2].x = D[2].x;
      B[2].y = D[2].y-WingRad;
      //      std::cout << "B2 " << Origin+X*B[2].x - Y*B[2].y << std::endl; 


      dirX = X;
      sideAngle[2] = atan2(C[1].x-A[2].x, C[1].y-A[2].y)*180/M_PI;
      Geometry::Quaternion::calcQRotDeg(sideAngle[2], Z).rotate(dirX);
      //      sideWingPlane[2] = ModelSupport::buildPlane(SMap, SI+8, Origin+Geometry::Vec3D(A[2].x, A[2].y, 0), dirX);
      sideWingPlane[2] = ModelSupport::buildPlane(SMap, SI+8, Origin+X*A[2].x -Y*A[2].y, dirX);

      // virtural plane
      dirX = X;
      virtAngle[2] = atan2(A[2].x-B[2].x, A[2].y-B[2].y)*180/M_PI;
      Geometry::Quaternion::calcQRotDeg(virtAngle[2], Z).rotate(dirX);
      //      ModelSupport::buildPlane(SMap, SI+30, Origin+Geometry::Vec3D(A[2].x, A[2].y, 0), dirX);
      ModelSupport::buildPlane(SMap, SI+30, Origin+X*A[2].x -Y*A[2].y, dirX);



      // cylinder in the x-y- corner
      D[3].x = D[0].x;
      D[3].y = D[2].y;
      ModelSupport::buildCylinder(SMap, SI+39, Origin+X*(D[3].x)-Y*(D[3].y), Z, WingRad);
      //      std::cout << "c39 " << Origin+X*(D[3].x)-Y*(D[3].y) << std::endl; // OK

      
      // intersect of cylinder 39 and plane 2
      A[3].x = A[0].x;
      A[3].y = A[2].y;
      //      std::cout << "A3 " << Origin+X*A[3].x - Y*A[3].y << std::endl; // ok
      
      // intersect of cylinder 39 and plane 3
      B[3].x = D[3].x;
      B[3].y = B[2].y;
      std::cout << "B3 " << Origin+X*B[3].x - Y*B[3].y << std::endl; 
      
      
      dirX = X;
      sideAngle[3] = atan2(C[1].x-A[3].x, C[1].y-A[3].y)*180/M_PI;
      Geometry::Quaternion::calcQRotDeg(sideAngle[3], Z).rotate(dirX);
      sideWingPlane[3] = ModelSupport::buildPlane(SMap, SI+2, Origin+X*A[3].x -Y*A[3].y, dirX);

      // virtural plane
      dirX = X;
      virtAngle[3] = atan2(A[3].x-B[3].x, A[3].y-B[3].y)*180/M_PI;
      Geometry::Quaternion::calcQRotDeg(virtAngle[3], Z).rotate(dirX);
      ModelSupport::buildPlane(SMap, SI+40, Origin+X*A[3].x - Y*A[3].y, dirX);

      // Roundings of the noses
      double theta3 = M_PI/2+atan((C[0].y-A[0].y)/(C[0].x-A[0].x)); // angle of inclination of the outer layer of plane 1
      D[4].x = C[0].x;
      D[4].y = C[0].y + WingCentralRad/sin(theta3);
      ModelSupport::buildCylinder(SMap, SI+46, Origin+X*(D[4].x)-Y*(D[4].y), Z, WingCentralRad);
      A[4].x = D[4].x - WingCentralRad*cos(theta3);
      A[4].y = D[4].y - WingCentralRad*sin(theta3);
      A[5].x = D[4].x + WingCentralRad*cos(theta3);
      A[5].y = A[4].y;
      ModelSupport::buildPlane(SMap, SI+47, Origin-Y*(A[4].y), Y);

      D[5].x = C[1].x;
      D[5].y = C[1].y - WingCentralRad/sin(theta3);
      ModelSupport::buildCylinder(SMap, SI+48, Origin+X*(D[5].x)-Y*(D[5].y), Z, WingCentralRad);
      A[6].x = D[5].x - WingCentralRad*cos(theta3);
      A[6].y = D[5].y + WingCentralRad*sin(theta3);
      A[7].x = D[5].x + WingCentralRad*cos(theta3);
      A[7].y = A[5].y;
      ModelSupport::buildPlane(SMap, SI+49, Origin-Y*(A[6].y), Y);

      SI += 50;
      //      fullWidth += Width[i];
    }
    WingsSeparated = D[4].y-WingCentralRad > D[5].y+WingCentralRad ? true : false;

    fullWidth = Width[nLayers-1];

    // Thermal moderator
    // E is the point where the outer layer of the cold moderator intersects the line y=0
    //    Geometry::Vec3D *ppNorm[4]; // normals

    ModelSupport::buildPlane(SMap, SI+5, Origin, X); // plane along the y-axis to join thermal layers x- and x+ when  ModWingLengh < (ModLength/2-total_layer_thickness)

    double h = (Height[nLayers-1]-Height[nLayers-2])/2.0; // thickness of the outer layer

    ModelSupport::buildPlane(SMap, SI+6, Origin-Z*(Height[nLayers-2]/2.0+PreBottomHeight), Z); // zmin of water
    ModelSupport::buildPlane(SMap, SI+7, Origin+Z*(Height[nLayers-2]/2.0+PreTopHeight+h), Z); // zmax of water
    ModelSupport::buildPlane(SMap, SI+8, Origin-Z*(Height[nLayers-2]/2.0+PreBottomHeight+0.0*PreWallThick), Z); // zmin of water + wall thick*0 since we do not need Al there !!! todo: remove this surface at all and use surf 6 instead
    ModelSupport::buildPlane(SMap, SI+9, Origin+Z*(Height[nLayers-2]/2.0+PreTopHeight+h+PreWallThick), Z); // zmax of water + wall thick

    Geometry::Plane *p;
    // x<0
    E[0].x = (C[0].x-A[0].x)*(0-A[0].y)/(C[0].y-A[0].y)+A[0].x;
    E[0].y = 0;

    // x>0
    E[1].x = ((C[1].x-A[2].x)*(0-A[2].y)/(C[1].y-A[2].y)+A[2].x);
    E[1].y = 0;

    double theta3 = M_PI/2+atan((C[0].y-A[0].y)/(C[0].x-A[0].x)); // angle of inclination of the outer layer of plane 1
    //    std::cout << "theta3 "  << theta3*180/M_PI << std::endl;
    
    const double xmin1 = (E[0].x < 0 ? E[0].x : 0.0); 
    const double xmin2 = (E[1].x > 0 ? E[1].x : 0.0);
    // x-y+
    // F is the point at the inclined planes
    F[0].x = xmin1 - PreLength*sin(theta3); // if the hydrogen halfs (x<0 and x>0) are split then start counting from zero
    F[0].y = (C[0].y<0 ? E[0].y : C[0].y) + PreLength*cos(theta3); // - / - / - then start counting Y from C[0].y

    G[0].x = F[0].x - PreThick1*cos(theta3); // yes, here is cos, and in Fx is sin - that's correct because EF is perp to GF
    G[0].y = F[0].y - PreThick1*sin(theta3);

    p = ModelSupport::buildPlane(SMap, SI+1, Origin+X*F[0].x-Y*F[0].y, Origin+X*G[0].x - Y*G[0].y, Origin+X*G[0].x - Y*G[0].y + Z*10, Origin+X*D[0].x -Y*D[0].y);
    ModelSupport::buildPlane(SMap, SI+2, Origin+X*(F[0].x-PreWallThick*sin(theta3)) -Y*(F[0].y+PreWallThick*cos(theta3)), p->getNormal());
    
    p = ModelSupport::buildPlane(SMap, SI+3, Origin+X*G[0].x-Y*G[0].y, Origin+X*(xmin1-PreBaseThick)-Y*E[0].y, Origin+X*G[0].x - Y*G[0].y + Z*10,
				 Origin-X*D[3].x - Y*D[3].y);
    ModelSupport::buildPlane(SMap, SI+4, Origin+X*(G[0].x-PreWallThick*sin(theta3)) - Y*(G[0].y-PreWallThick*cos(theta3)), p->getNormal());

    // water undeneath the cold part
    K[0].x = F[0].x+PreThick2*cos(theta3);
    K[0].y = F[0].y+PreThick2*sin(theta3);

    ModelSupport::buildPlane(SMap, SI+53, Origin+X*K[0].x - Y*K[0].y, sideWingPlane[0]->getNormal()); // water
    ModelSupport::buildPlane(SMap, SI+54, Origin+X*(K[0].x+PreWallThick*cos(theta3)) - Y*(K[0].y+PreWallThick*cos(theta3)), sideWingPlane[0]->getNormal()); // wall

    // x-y-
    theta3 = atan((C[1].x-A[3].x)/(C[1].y-A[3].y));
    //    std::cout << "theta3 "  << theta3*180/M_PI << std::endl;
    F[3].x = xmin1 - PreLength*sin(theta3);
    F[3].y = (C[1].y>0 ? E[0].y : C[1].y) - PreLength*cos(theta3);
    // std::cout << "F0: " << F[0].x << " " << F[0].y << std::endl;
    // std::cout << "G0: " << G[0].x << " " << G[0].y << std::endl;
    // std::cout << "F3: " << F[3].x << " " << F[3].y << std::endl;

    G[3].x = F[3].x - PreThick1*cos(theta3);
    G[3].y = F[3].y + PreThick1*sin(theta3);

    p = ModelSupport::buildPlane(SMap, SI+11, Origin+X*F[3].x - Y*F[3].y + Z*0, Origin+X*G[3].x - Y*G[3].y + Z*0, Origin+X*G[3].x - Y*G[3].y + Z*10, Origin+X*D[3].x - Y*D[3].y + Z*0);
    ModelSupport::buildPlane(SMap, SI+12, Origin+X*(G[3].x-PreWallThick*sin(theta3)) - Y*(G[3].y-PreWallThick*cos(theta3)), p->getNormal());

    p = ModelSupport::buildPlane(SMap, SI+13, Origin+X*(G[3].x) - Y*(G[3].y), Origin+X*(xmin1-PreBaseThick) - Y*(E[0].y), Origin+X*(G[3].x) - Y*(G[3].y) + Z*(10), Origin+X*(D[0].x) - Y*(D[0].y));
    ModelSupport::buildPlane(SMap, SI+14, Origin+X*(G[3].x-PreWallThick*sin(theta3)) - Y*(G[3].y+PreWallThick*cos(theta3)), p->getNormal());

    // water undeneath the cold part
    K[3].x = F[3].x+PreThick2*cos(theta3);
    K[3].y = F[3].y-PreThick2*sin(theta3);
    ModelSupport::buildPlane(SMap, SI+83, Origin+X*K[3].x - Y*K[3].y, sideWingPlane[3]->getNormal()); // water
    ModelSupport::buildPlane(SMap, SI+84, Origin+X*(K[3].x+PreWallThick*cos(theta3)) - Y*(K[3].y-PreWallThick*cos(theta3)), sideWingPlane[3]->getNormal()); // wall


    // x+y-
    theta3 = -atan((C[1].x-A[2].x)/(C[1].y-A[2].y));
    F[2].x = xmin2 + PreLength*sin(theta3); // if the hydrogen halves are split, start counting from zero
    F[2].y = (C[1].y > 0 ? E[1].y : C[1].y) - PreLength*cos(theta3);
    //    std::cout << "F2: " << F[2].x << " " << F[2].y << " " << theta3 << std::endl;

    G[2].x = F[2].x + PreThick1*cos(theta3);
    G[2].y = F[2].y + PreThick1*sin(theta3);
    //    std::cout << "G2: " << G[2].x << " " << G[2].y << " " << theta3 << std::endl;
 
    p = ModelSupport::buildPlane(SMap, SI+21, Origin+X*(F[2].x) - Y*(F[2].y), Origin+X*(G[2].x) - Y*(G[2].y), Origin+X*(G[2].x) - Y*(G[2].y) + Z*(10), Origin+X*(D[2].x) - Y*(D[2].y)); // A1
    ModelSupport::buildPlane(SMap, SI+22, Origin+X*(G[2].x+PreWallThick*sin(theta3)) - Y*(G[2].y-PreWallThick*cos(theta3)), p->getNormal());

    p = ModelSupport::buildPlane(SMap, SI+23, Origin+X*(G[2].x) - Y*(G[2].y), Origin+X*(xmin2+PreBaseThick) - Y*(E[1].y), Origin+X*(G[2].x) - Y*(G[2].y) + Z*(10), Origin+X*(A[1].x) - Y*(A[1].y));
    ModelSupport::buildPlane(SMap, SI+24, Origin+X*(G[2].x+PreWallThick*sin(theta3)) - Y*(G[2].y+PreWallThick*cos(theta3)), p->getNormal());

    // water undeneath the cold part
    K[2].x = F[2].x-PreThick2*cos(theta3);
    K[2].y = F[2].y-PreThick2*sin(theta3);
    ModelSupport::buildPlane(SMap, SI+93, Origin+X*K[2].x - Y*K[2].y, sideWingPlane[2]->getNormal()); // water
    ModelSupport::buildPlane(SMap, SI+94, Origin+X*(K[2].x-PreWallThick*cos(theta3)) - Y*(K[2].y-PreWallThick*cos(theta3)), sideWingPlane[2]->getNormal()); // wall


    // x+y+
    theta3 = atan((C[0].x-A[1].x)/(C[0].y-A[1].y));
    F[1].x = F[2].x; //xmin2 + PreLength*sin(theta3);
    F[1].y = (C[0].y<0 ? E[1].y : C[0].y) + PreLength*cos(theta3);

    G[1].x = F[1].x + PreThick1*cos(theta3);
    G[1].y = F[1].y - PreThick1*sin(theta3);
    //    std::cout << "G1: " << G[1].x << " " << G[1].y << " " << theta3 << std::endl;


    p = ModelSupport::buildPlane(SMap, SI+31, Origin+X*(F[1].x) - Y*(F[1].y), Origin+X*(G[1].x) - Y*(G[1].y), Origin+X*(G[1].x) - Y*(G[1].y) + Z*(10), Origin+X*(D[1].x) - Y*(D[1].y)); // A2

    //    cylph << "F5xb " << F[1].x+PreWallThick*sin(theta3) << "\nF5yb " << F[1].y+PreWallThick*cos(theta3) << std::endl;

    ModelSupport::buildPlane(SMap, SI+32, Origin+X*(G[1].x+PreWallThick*sin(theta3)) - Y*(G[1].y+PreWallThick*cos(theta3)), p->getNormal());

    p = ModelSupport::buildPlane(SMap, SI+33, Origin+X*(G[1].x) - Y*(G[1].y), Origin+X*(xmin2+PreBaseThick) - Y*(E[1].y), Origin+X*(G[1].x) - Y*(G[1].y) + Z*(10), Origin+X*(A[2].x) - Y*(A[2].y));
    ModelSupport::buildPlane(SMap, SI+34, Origin+X*(G[1].x+PreWallThick*sin(theta3)) - Y*(G[1].y-PreWallThick*cos(theta3)), p->getNormal());

    // water undeneath the cold part
    K[1].x = F[1].x-PreThick2*cos(theta3);
    K[1].y = F[1].y+PreThick2*sin(theta3);
    ModelSupport::buildPlane(SMap, SI+1103, Origin+X*K[1].x - Y*K[1].y, sideWingPlane[1]->getNormal()); // water
    ModelSupport::buildPlane(SMap, SI+1104, Origin+X*(K[1].x-PreWallThick*cos(theta3))-Y*(K[1].y+PreWallThick*cos(theta3)), sideWingPlane[1]->getNormal()); // wall

   
    // top pre moderator
    ModelSupport::buildPlane(SMap, SI+55, Origin+Z*(Height[nLayers-1]/2.0+TopPreHeight), Z); // see also plane 75
    if (TopPreType==1) {
      ModelSupport::buildPlane(SMap, SI+56, Origin-X*(Width[nLayers-2]/2), X);
      ModelSupport::buildPlane(SMap, SI+57, Origin+X*(Width[nLayers-2]/2), X);
      // virtual planes (like 10,20,30,40):
      dirX = X; Geometry::Quaternion::calcQRotDeg(45, Z).rotate(dirX);
      ModelSupport::buildPlane(SMap, SI+58, Origin+X*(D[0].x) - Y*(D[0].y+vWingRad[nLayers-2]), dirX); // as 10
      ModelSupport::buildPlane(SMap, SI+59, Origin+X*(D[2].x) - Y*(D[2].y-vWingRad[nLayers-2]), dirX); // as 30

      ModelSupport::buildPlane(SMap, SI+78, Origin+X*(D[0].x) - Y*(D[0].y+vWingRad[nLayers-1]), dirX); // as 10 + Al
      ModelSupport::buildPlane(SMap, SI+79, Origin+X*(D[2].x) - Y*(D[2].y-vWingRad[nLayers-1]), dirX); // as 30 + Al
      dirX = X; Geometry::Quaternion::calcQRotDeg(-45, Z).rotate(dirX);
      ModelSupport::buildPlane(SMap, SI+60, Origin+X*(D[1].x) - Y*(D[1].y+vWingRad[nLayers-2]), dirX); // as 20
      ModelSupport::buildPlane(SMap, SI+61, Origin+X*(D[3].x) - Y*(D[3].y-vWingRad[nLayers-2]), dirX); // as 40
      ModelSupport::buildPlane(SMap, SI+80, Origin+X*(D[1].x) - Y*(D[1].y+vWingRad[nLayers-1]), dirX); // as 20 + Al
      ModelSupport::buildPlane(SMap, SI+81, Origin+X*(D[3].x) - Y*(D[3].y-vWingRad[nLayers-1]), dirX); // as 40 + Al
    } else if (TopPreType==2) {
      ModelSupport::buildPlane(SMap, SI+56, Origin-X*(TopPreWidth/2), X);
      ModelSupport::buildPlane(SMap, SI+57, Origin+X*(TopPreWidth/2), X);
    } else if (TopPreType==3) { // cylindrical
      ModelSupport::buildCylinder(SMap, SI+56, Origin, Z, TopPreRadius);
      ModelSupport::buildCylinder(SMap, SI+57, Origin, Z, TopPreRadius+PreWallThick);

      ModelSupport::buildCylinder(SMap, SI+58, Origin, Z, TopPreRadius+PreWallThick+TopPreTopVoidHeight);
      ModelSupport::buildCylinder(SMap, SI+59, Origin, Z, TopPreRadius+PreWallThick+TopPreTopVoidHeight+TopPreTopAlHeight);
    }

    ModelSupport::buildPlane(SMap, SI+65, Origin+Z*(Height[nLayers-1]/2.0+TopPreHeight+PreWallThick), Z);
    if (TopPreType==1) {
      ModelSupport::buildPlane(SMap, SI+66, Origin-X*(Width[nLayers-1]/2), X);
      ModelSupport::buildPlane(SMap, SI+67, Origin+X*(Width[nLayers-1]/2), X);
    } else if (TopPreType==2) {
      ModelSupport::buildPlane(SMap, SI+66, Origin-X*(TopPreWidth/2+PreWallThick), X);
      ModelSupport::buildPlane(SMap, SI+67, Origin+X*(TopPreWidth/2+PreWallThick), X);
    }

    ModelSupport::buildPlane(SMap, SI+75, Origin+Z*(Height[nLayers-1]/2.0-PreWallThick), Z); // see also plane 55

    // second Al layer on top of TopPre
    ModelSupport::buildPlane(SMap, SI+76, Origin+Z*(Height[nLayers-1]/2.0+TopPreHeight+PreWallThick+TopPreTopVoidHeight), Z);
    ModelSupport::buildPlane(SMap, SI+176, Origin+Z*(Height[nLayers-1]/2.0+TopPreHeight+PreWallThick+TopPreTopVoidHeight+TopPreTopAlHeight), Z);




    // Flight lines
    /// x<0
    Geometry::Vec3D yDirc(Y);
    Geometry::Quaternion::calcQRotDeg(AFlightAngleXY1, Z).rotate(yDirc);
    ModelSupport::buildPlane(SMap, SI+103, Origin-Y*(Length[nLayers-1]/2.0) - X*(fullWidth/2.0-FlightXOffset), yDirc);
    ModelSupport::buildPlane(SMap, SI+113, Origin-Y*(Length[nLayers-1]/2.0) - X*(fullWidth/2.0-FlightXOffset - FlightLineWallThick/sin(AFlightAngleXY1*M_PI/180)), yDirc);
    if (FlightLineType==1) {
      ModelSupport::buildPlane(SMap, SI+143, Origin-Y*(Length[nLayers-1]/2.0) - X*(fullWidth/2.0-FlightXOffset), X);
      ModelSupport::buildPlane(SMap, SI+153, Origin-Y*(Length[nLayers-1]/2.0) - X*(fullWidth/2.0-FlightXOffset - FlightLineWallThick/sin(AFlightAngleXY1*M_PI/180)), X);

      ModelSupport::buildPlane(SMap, SI+144, Origin+Y*(Length[nLayers-1]/2.0) + X*(fullWidth/2.0-FlightXOffset), X);
      ModelSupport::buildPlane(SMap, SI+154, Origin+Y*(Length[nLayers-1]/2.0) + X*(fullWidth/2.0-FlightXOffset - FlightLineWallThick/sin(AFlightAngleXY1*M_PI/180)), X);
    }

    yDirc = Y;
    Geometry::Quaternion::calcQRotDeg(-AFlightAngleXY2, Z).rotate(yDirc);
    ModelSupport::buildPlane(SMap, SI+104, Origin+Y*(Length[nLayers-1]/2.0) - X*(fullWidth/2.0-FlightXOffset), yDirc);
    ModelSupport::buildPlane(SMap, SI+114, Origin+Y*(Length[nLayers-1]/2.0) - X*(fullWidth/2.0-FlightXOffset - FlightLineWallThick/sin(AFlightAngleXY2*M_PI/180)), yDirc);

    Geometry::Vec3D zDirc(Z);
    Geometry::Quaternion::calcQRotDeg(-AFlightAngleZBase, Y).rotate(zDirc);
    ModelSupport::buildPlane(SMap, SI+105, Origin-Z*(Height[nLayers-2]/2.0), zDirc); // here
    ModelSupport::buildPlane(SMap, SI+115, Origin-Z*(Height[nLayers-2]/2.0 + FlightLineWallThick), zDirc);
    zDirc = Z;
    Geometry::Quaternion::calcQRotDeg(AFlightAngleZTop, Y).rotate(zDirc);

    double zshift = 0.0;
    if (FlightLineWrapTopPreType==0)
      zshift = TopPreHeight + PreWallThick;
    if (TopPreType) {
      ModelSupport::buildPlane(SMap, SI+106, Origin+Z*(Height[nLayers-2]/2.0 + zshift), zDirc);
      ModelSupport::buildPlane(SMap, SI+116, Origin+Z*(Height[nLayers-2]/2.0 + FlightLineWallThick + zshift), zDirc);
      ModelSupport::buildPlane(SMap, SI+216, Origin+Z*(Height[nLayers-2]/2.0 + FlightLineWallThick + zshift+TopPreTopVoidHeight), zDirc);
      ModelSupport::buildPlane(SMap, SI+316, Origin+Z*(Height[nLayers-2]/2.0 + FlightLineWallThick + zshift+TopPreTopVoidHeight+TopPreTopAlHeight), zDirc);
      } else {
      ModelSupport::buildPlane(SMap, SI+106, Origin+Z*(Height[0]/2.0), zDirc);
      ModelSupport::buildPlane(SMap, SI+116, Origin+Z*(Height[0]/2.0+FlightLineWallThick), zDirc);
    }
    // x>0
    yDirc = Y;
    Geometry::Quaternion::calcQRotDeg(-BFlightAngleXY1, Z).rotate(yDirc);
    //    ModelSupport::buildPlane(SMap, SI+123, Origin-Y*(Length[nLayers-1]/2.0) + X*(fullWidth/2.0-WingRad), yDirc);
    ModelSupport::buildPlane(SMap, SI+123, Origin-Y*(Length[nLayers-1]/2.0) + X*(fullWidth/2.0-FlightXOffset), yDirc);
    //    ModelSupport::buildPlane(SMap, SI+133, Origin-Y*(Length[nLayers-1]/2.0) + X*(fullWidth/2.0-WingRad  -FlightLineWallThick/sin(BFlightAngleXY1*M_PI/180)), yDirc);
    ModelSupport::buildPlane(SMap, SI+133, Origin-Y*(Length[nLayers-1]/2.0) + X*(fullWidth/2.0-FlightXOffset  -FlightLineWallThick/sin(BFlightAngleXY1*M_PI/180)), yDirc);
    yDirc = Y;
    Geometry::Quaternion::calcQRotDeg(BFlightAngleXY2, Z).rotate(yDirc);
    //    ModelSupport::buildPlane(SMap, SI+124, Origin+Y*(Length[nLayers-1]/2.0) + X*(fullWidth/2.0-WingRad), yDirc);
    ModelSupport::buildPlane(SMap, SI+124, Origin+Y*(Length[nLayers-1]/2.0) + X*(fullWidth/2.0-FlightXOffset), yDirc);
    //    ModelSupport::buildPlane(SMap, SI+134, Origin+Y*(Length[nLayers-1]/2.0) + X*(fullWidth/2.0-WingRad - FlightLineWallThick/sin(BFlightAngleXY2*M_PI/180.0)), yDirc);
    ModelSupport::buildPlane(SMap, SI+134, Origin+Y*(Length[nLayers-1]/2.0) + X*(fullWidth/2.0-FlightXOffset - FlightLineWallThick/sin(BFlightAngleXY2*M_PI/180.0)), yDirc);

    zDirc = Z;
    Geometry::Quaternion::calcQRotDeg(BFlightAngleZBase, Y).rotate(zDirc);
    ModelSupport::buildPlane(SMap, SI+125, Origin-Z*(Height[nLayers-2]/2.0), zDirc); // here
    ModelSupport::buildPlane(SMap, SI+135, Origin-Z*(Height[nLayers-2]/2.0 + FlightLineWallThick), zDirc);
    zDirc = Z;
    Geometry::Quaternion::calcQRotDeg(-BFlightAngleZTop, Y).rotate(zDirc);

    if (TopPreType) {
      ModelSupport::buildPlane(SMap, SI+126, Origin+Z*(Height[nLayers-2]/2.0 + zshift), zDirc);
      ModelSupport::buildPlane(SMap, SI+136, Origin+Z*(Height[nLayers-2]/2.0 + FlightLineWallThick + zshift), zDirc);
    } else {
      ModelSupport::buildPlane(SMap, SI+126, Origin+Z*(Height[0]/2.0), zDirc);
      ModelSupport::buildPlane(SMap, SI+136, Origin+Z*(Height[0]/2.0 + FlightLineWallThick), zDirc);
    }

    //    cylph.close();
    return; 
  }

  void
  ButterflyModerator::addToInsertChain(attachSystem::ContainedComp& CC) const
  /*!
    Adds this object to the containedComp to be inserted.
    \param CC :: ContainedComp object to add to this
  */
  {
    for(int i=modIndex+1;i<cellIndex;i++) CC.addInsertCell(i);
    
    return;
  }


  void
  ButterflyModerator::createObjects(Simulation& System,
				    const attachSystem::FixedComp& ShutterBay, const long int sIndex)
  {
    // Builds the moderator
    // see archetypes/ButterflyModerator

    ELog::RegMethod RegA("ButterflyModerator","createObjects");

    // side surface
    const std::string sSurf=(sIndex<0) ? ShutterBay.getLinkComplement(static_cast<size_t>(-(sIndex+1))) : ShutterBay.getLinkString(static_cast<size_t>(sIndex));

    std::string Out, Out1, Out2;
    HeadRule HR;
    int SI(modIndex);

    for (size_t i=0; i<nLayers; i++) {
      //                                                 right wing                   left wing             round corners
      // no nose rounding:     
      if (WingsSeparated) // round the nose
	Out = ModelSupport::getComposite(SMap, SI, " (((1 3 -11 -7 10 -20 -47) : (2 -4 -8 -30 40 11 49) : -9 : -19 : -29 : -39 : -46 : -48) 5 -6) "); 
      else
	Out = ModelSupport::getComposite(SMap, SI, " (((1 3 -11 -7 10 -20     ) : (2 -4 -8 -30 40 11   ) : -9 : -19 : -29 : -39) 5 -6 ) "); 
      Out1 = Out;

      if (i) {
	//	Out += ModelSupport::getComposite(SMap, SI-50, " ((-1:-3:11:7:-10:20:-5:6) (-2:4:8:30:-40:-11:-5:6) (9:-5:6) (19:-5:6) (29:-5:6) (39:-5:6) ) "); 
	// next lines are the same as above with making use of HR
	//	HeadRule HR;
	HR.procString(Out2); // at this point Out2 = OUt1 at SI-50 offset
	HR.makeComplement();
	Out += HR.display();
      }

      System.addCell(MonteCarlo::Qhull(cellIndex++, mat[i], temp[i], Out));
      SI += 50;
      Out2 = Out1;
    }

    addOuterUnionSurf(Out1);

    // Thermal moderator
      if (WingsSeparated) { // round the nose
	Out = ModelSupport::getComposite(SMap, SI, SI-50, " -1 -53 -3 6 -7 -5 -11M 9M (-1M:-5M:6M) : (1M 46M 47M -5 6 -7 48M) ");
	//System.addCell(MonteCarlo::Qhull(cellIndex++, PreMat, PreTemp, Out));
	Out += ":" + ModelSupport::getComposite(SMap, SI, SI-50, " -11 -83 -13 6 -7 -5 11M 39M (-2M:-5M:6M) : (2M 48M -49M -5 6 -7 46M) "); 

	Out += ":" + ModelSupport::getComposite(SMap, SI, SI-50, " -21 93 -23 6 -7 5 11M 29M (8M:-5M:6M) : (-8M 48M -49M 5 6 -7 46M) ");
	//	System.addCell(MonteCarlo::Qhull(cellIndex++, PreMat, PreTemp, Out));
	Out += ":" + ModelSupport::getComposite(SMap, SI, SI-50, " -31 1103 -33 6 -7 5 -11M 19M (7M:-5M:6M) : (-7M 46M 47M 5 6 -7 48M) ");

	if (TopPreType) {replaceSurface(Out, SI, SI-50, "-7", "-6"); replaceSurface(Out, SI, SI-50, "7", "6");}
	System.addCell(MonteCarlo::Qhull(cellIndex++, PreMat, PreTemp, Out));
      }
      else {
	Out = ModelSupport::getComposite(SMap, SI, SI-50, "(-1 -53 -3 6 -7 -5 -11M 9M (-1M:-5M:6M)                       ) : (-11 -83 -13 6 -7 -5 11M 39M (-2M:-5M:6M)) ");
	Out += ":" + ModelSupport::getComposite(SMap, SI, SI-50, "(-21 93 -23 6 -7 5 11M 29M (8M:-5M:6M) ) : (-31 1103 -33 6 -7 5 -11M 19M (7M:-5M:6M) )");
	if (TopPreType) {replaceSurface(Out, SI, SI-50, "-7", "-6"); replaceSurface(Out, SI, SI-50, "7", "6");}
	System.addCell(MonteCarlo::Qhull(cellIndex++, PreMat, PreTemp, Out));
      }
      //      Out = ModelSupport::getComposite(SMap, SI, SI-50, "(-2 -4 -54 -11M -5 8 -9 9M (1:53:3:-6:7) (-1M:-5M:6M)) : (-12 -14 -84 11M -5 8 -9 39M (11:83:13:-6:7) (-2M:-5M:6M))");
      Out = ModelSupport::getComposite(SMap, SI, SI-50, "-2 -4 -54 -11M -5 8 -9 9M (1:53:3:-6:7) (-1M:-5M:6M)");
      if (TopPreType) {replaceSurface(Out, SI, SI-50, "-9", "-6"); replaceSurface(Out, SI, SI-50, "7", "6");}
      System.addCell(MonteCarlo::Qhull(cellIndex++, PreWallMat, PreWallTemp, Out));
      Out = ModelSupport::getComposite(SMap, SI, SI-50, "-12 -14 -84 11M -5 8 -9 39M (11:83:13:-6:7) (-2M:-5M:6M)");
      if (TopPreType) {replaceSurface(Out, SI, SI-50, "-9", "-6"); replaceSurface(Out, SI, SI-50, "7", "6");}
      System.addCell(MonteCarlo::Qhull(cellIndex++, PreWallMat, PreWallTemp, Out));
      Out = ModelSupport::getComposite(SMap, SI, SI-50, "(-2 -4 -54 -11M -5 8 -9 9M) : (-12 -14 -84 11M -5 8 -9 39M)");
      if (TopPreType) {replaceSurface(Out, SI, SI-50, "-9", "-6"); replaceSurface(Out, SI, SI-50, "7", "6");}
      addOuterUnionSurf(Out); 
      
    

      // water
      /*      if (WingsSeparated) {// round the nose
	Out = ModelSupport::getComposite(SMap, SI, SI-50, " -21 93 -23 6 -7 5 11M 29M (8M:-5M:6M) : (-8M 48M -49M 5 6 -7 46M) ");
	//	System.addCell(MonteCarlo::Qhull(cellIndex++, PreMat, PreTemp, Out));
	Out = Out + ":" + ModelSupport::getComposite(SMap, SI, SI-50, " -31 1103 -33 6 -7 5 -11M 19M (7M:-5M:6M) : (-7M 46M 47M 5 6 -7 48M) ");
	System.addCell(MonteCarlo::Qhull(cellIndex++, PreMat, PreTemp, Out));
      }  else {
	Out = ModelSupport::getComposite(SMap, SI, SI-50, "(-21 93 -23 6 -7 5 11M 29M (8M:-5M:6M) ) : (-31 1103 -33 6 -7 5 -11M 19M (7M:-5M:6M) )");
	System.addCell(MonteCarlo::Qhull(cellIndex++, PreMat, PreTemp, Out));
	}*/
      //      addOuterUnionSurf(Out); 
      // Al
      Out = ModelSupport::getComposite(SMap, SI, SI-50, "-22 -24 94 11M 5 8 -9 29M (21:-93:23:-6:7) (8M:-5M:6M)");
      if (TopPreType) {replaceSurface(Out, SI, SI-50, "-9", "-6"); replaceSurface(Out, SI, SI-50, "7", "6");}
      System.addCell(MonteCarlo::Qhull(cellIndex++, PreWallMat, PreWallTemp, Out));
      Out = ModelSupport::getComposite(SMap, SI, SI-50, "-32 -34 1104 -11M 5 8 -9 19M (31:-1103:33:-6:7) (7M:-5M:6M)");
      if (TopPreType) {replaceSurface(Out, SI, SI-50, "-9", "-6"); replaceSurface(Out, SI, SI-50, "7", "6");}
      System.addCell(MonteCarlo::Qhull(cellIndex++, PreWallMat, PreWallTemp, Out));
      Out = ModelSupport::getComposite(SMap, SI, SI-50, "(-22 -24 94 11M 5 8 -9 29M) : (-32 -34 1104 -11M 5 8 -9 19M)");
      if (TopPreType) {replaceSurface(Out, SI, SI-50, "-9", "-6");}
      addOuterUnionSurf(Out); 

      // complementary mod surface
      // !!! todo: to make tracing faster, substitute from the flight lines only the relevant part of the moderator, but not its full volume
      //     but from another point of view substituting all volume is safer since it allows any arbitrary geometry of flight lines
      std::string ColdWingsComplement = ModelSupport::getComposite(SMap, SI-50, " (-1:-3:11:7:-10:20:-5:6)  (-2:4:8:30:-40:-11:-5:6) (9:-5:6) (19:-5:6) (29:-5:6) (39:-5:6) ");
      //                                                                                       water vessel at x<0                            water vessel at x>0
      std::string WaterWingsComplement = ModelSupport::getComposite(SMap, SI, SI-50, "(2:4:1M:11M:-8:9:5) (12:14:2M:-11M:-8:9:5) (22:24:-8M:-11M:-8:9:-5) (32:34:-7M:11M:-8:9:-5)");
      std::string ModComplement = ColdWingsComplement + WaterWingsComplement;


      
      // TopPre
      std::string TopPreOuterSurface("");
      if (TopPreType) {
	// water
	Out = ModelSupport::getComposite(SMap, SI, SI-50, " -55 6M ");
	if ((TopPreType==1) or (TopPreType==2))
	  Out += ModelSupport::getComposite(SMap, SI, SI-50, " (56 -57  ");
	if (TopPreType==1)
	  Out += ModelSupport::getComposite(SMap, SI-100, SI, " 58M -59M -60M 61M 3 -4 : -9 : -19 : -29 : -39)");
	else if (TopPreType==2)
	  Out += ReflectorSideAl + ")";
	if (TopPreType==3) {
	  Out += ModelSupport::getComposite(SMap, SI, " -56 ");
	}
	System.addCell(MonteCarlo::Qhull(cellIndex++, PreMat, PreTemp, Out));

	// Al
	Out1 = Out;
	HR.procString(Out1);
	HR.makeComplement();

	Out = ModelSupport::getComposite(SMap, SI, SI-50, " -65 75 ");
	if ((TopPreType==1) || (TopPreType==2)) {
	  Out += ModelSupport::getComposite(SMap, SI, SI-50, " (66 -67 ");
	  if (TopPreType==1) 
	    Out += ModelSupport::getComposite(SMap, SI-50, SI, " 78M -79M -80M 81M 3 -4 : -9 : -19 : -29 : -39)"); // 
	  else if (TopPreType==2)
	    Out += " +" + ReflectorSideBe + ")"; 
	} else if (TopPreType==3) {
	  Out += ModelSupport::getComposite(SMap, SI, SI-50, " -65 75 -57 ");
	}
	addOuterUnionSurf(Out); 
	Out1 = WaterWingsComplement;
	replaceSurface(Out1, SI, SI-50, "9", "6M"); // replace upper surface of the cross by the upper surface of the cold wing
	replaceSurface(Out1, SI, SI, "8", "75"); // just to gain speed: replace lower water surface by the lower surface of TopPreMod al container
	Out += ColdWingsComplement + Out1;
	TopPreOuterSurface = Out;
	Out += HR.display();
	System.addCell(MonteCarlo::Qhull(cellIndex++, PreWallMat, PreWallTemp, Out));
	setCell(keyName+"Wall",1,cellIndex-1); // for TSupply Pipe - name this cell in order to remove it in makeESS: TopSupplyPipe->addInsertCell
	//	ELog::EM << "SET CELL : " << keyName << "Wall " << cellIndex-1 << " " << TopPreType << ELog::endCrit;

	//setCell(keyName+"Ring",1,cellIndex-1); // for TSupply Pipe - name this cell in order to remove it in makeESS: TopSupplyPipe->addInsertCell

	// void+Al above Butterfly and TopPre
	if (1) {
	  if (TopPreType==1) { // todo: remove (1) and add in Table
	    Out = ModelSupport::getComposite(SMap, SI, SI-50, " 65 -76 (66 -67 78 -79 -80 81 3M -4M : -9M : -19M : -29M : -39M) ");
	    System.addCell(MonteCarlo::Qhull(cellIndex++, 0, 0, Out));
	    addOuterUnionSurf(Out);

	    Out = ModelSupport::getComposite(SMap, SI, SI-50, " 76 -176 (66 -67 78 -79 -80 81 3M -4M : -9M : -19M : -29M : -39M) ");
	    System.addCell(MonteCarlo::Qhull(cellIndex++, PreWallMat, PreWallTemp, Out));
	    addOuterUnionSurf(Out);
	  } else if (TopPreType==3) {
	    //	    ELog::EM << "here" << ELog::endErr;
	    Out = ModelSupport::getComposite(SMap, SI, SI-50, " 65 -76 -58 "); // void aboive TopPre
	    System.addCell(MonteCarlo::Qhull(cellIndex++, 0, 0, Out));
	    addOuterUnionSurf(Out);

	    Out = ModelSupport::getComposite(SMap, SI, SI-50, " -65 116 -58 57 "); // void outside rad TopPre
	    System.addCell(MonteCarlo::Qhull(cellIndex++, 0, 0, Out));
	    addOuterUnionSurf(Out);

	    Out = ModelSupport::getComposite(SMap, SI, SI-50, " 76 -176 -59 "); // Al above this void
	    System.addCell(MonteCarlo::Qhull(cellIndex++, PreWallMat, PreWallTemp, Out));
	    addOuterUnionSurf(Out);

	    Out = ModelSupport::getComposite(SMap, SI, SI-50, " -76 216 -59 58 "); // Al above this void
	    System.addCell(MonteCarlo::Qhull(cellIndex++, PreWallMat, PreWallTemp, Out));
	    addOuterUnionSurf(Out);
	  }
	}
	
	HR.procString(TopPreOuterSurface);
	HR.makeComplement();
	TopPreOuterSurface = HR.display();
      }
      ModComplement += TopPreOuterSurface;
      if (TopPreType) replaceSurface(ModComplement, SI, SI-50, "9", "6");
      
    // exclude thermal moderator
    //    ModComplement += ModelSupport::getComposite(SMap, SI, SI-50, "(-1:1M:11M:-5M:6M)");

      // water in BeRef at the level of Butterfly:
      Out = ModelSupport::getComposite(SMap, SI, SI-50, " 105 -106 -3M ") + ReflectorSideAl;
      if (FlightLineType==1) Out += ModelSupport::getComposite(SMap, SI, " 153 -154 ");
      else Out += ModelSupport::getComposite(SMap, SI, " -113 -133 ");
      System.addCell(MonteCarlo::Qhull(cellIndex++, MaterialInBeRef, 0.0, Out)); 
      addOuterUnionSurf(Out);

      Out = ModelSupport::getComposite(SMap, SI, SI-50, " 105 -106 4M ") + ReflectorSideAl;
      if (FlightLineType==1) Out += ModelSupport::getComposite(SMap, SI, " 153 -154 ");
      else Out += ModelSupport::getComposite(SMap, SI, " 114 134 ");
      System.addCell(MonteCarlo::Qhull(cellIndex++, MaterialInBeRef, 0.0, Out)); 
      addOuterUnionSurf(Out);
    
      // top Al + void + Al
      if (1) {
	Out1 = ReflectorSideAl;
	if (TopPreType==3) Out1 += ModelSupport::getComposite(SMap, SI, " 57 ");
	if ((TopPreType==1) || (TopPreType==3)) { // todo: uncomment (1) after added in the Table. optimise addOuterUnionSurf
	  // Al
	  Out = ModelSupport::getComposite(SMap, SI, SI-50, " 106 -116 -3M ") + Out1; // -3M
	  if (FlightLineType==1) Out += ModelSupport::getComposite(SMap, SI, " 153 -154 ");
	  else Out += ModelSupport::getComposite(SMap, SI, " -113 -133 ");
	  System.addCell(MonteCarlo::Qhull(cellIndex++, PreWallMat, PreWallTemp, Out)); 

	  Out = ModelSupport::getComposite(SMap, SI, SI-50, " 106 -116 4M ") + Out1; // 4M
	  if (FlightLineType==1) Out += ModelSupport::getComposite(SMap, SI, " 153 -154 ");
	  else Out += ModelSupport::getComposite(SMap, SI, " 114 134 ");
	  System.addCell(MonteCarlo::Qhull(cellIndex++, PreWallMat, PreWallTemp, Out)); 

	  // void layer
	  Out = ModelSupport::getComposite(SMap, SI, SI-50, " 116 -216 58 ") + Out1; // -3M 
	  //	  if (FlightLineType==1) Out += ModelSupport::getComposite(SMap, SI, " 153 -154 ");
	  //	  else Out += ModelSupport::getComposite(SMap, SI, " -113 -133 ");
	  System.addCell(MonteCarlo::Qhull(cellIndex++, 0, 0.0, Out)); 

	  //	  Out = ModelSupport::getComposite(SMap, SI, SI-50, " 116 -216 ") + Out1; // 4M 
	  //	  if (FlightLineType==1) Out += ModelSupport::getComposite(SMap, SI, " 153 -154 ");
	  //	  else Out += ModelSupport::getComposite(SMap, SI, " 114 134 ");
	  //	  System.addCell(MonteCarlo::Qhull(cellIndex++, 0, 0.0, Out)); 

	  // Al
	  Out = ModelSupport::getComposite(SMap, SI, SI-50, " 216 -316 59 ") + Out1; // -3M 
	  //	  if (FlightLineType==1) Out += ModelSupport::getComposite(SMap, SI, " 153 -154 ");
	  //	  else Out += ModelSupport::getComposite(SMap, SI, " -113 -133 ");
	  System.addCell(MonteCarlo::Qhull(cellIndex++, PreWallMat, PreWallTemp, Out)); 

	  //	  Out = ModelSupport::getComposite(SMap, SI, SI-50, " 216 -316 ") + Out1; // 4M 
	  //	  if (FlightLineType==1) Out += ModelSupport::getComposite(SMap, SI, " 153 -154 ");
	  //	  else Out += ModelSupport::getComposite(SMap, SI, " 114 134 ");
	  //	  System.addCell(MonteCarlo::Qhull(cellIndex++, PreWallMat, PreWallTemp, Out)); 


	  // outer surface
	  Out = ModelSupport::getComposite(SMap, SI, SI-50, " 106 -316 ") + Out1; // -3M 
	  //	  if (FlightLineType==1) Out += ModelSupport::getComposite(SMap, SI, " 153 -154 ");
	  //	  else Out += ModelSupport::getComposite(SMap, SI, " -113 -133 ");
	  addOuterUnionSurf(Out);

	  Out = ModelSupport::getComposite(SMap, SI, SI-50, " 106 -316 ") + Out1; // 4M 
	  //	  if (FlightLineType==1) Out += ModelSupport::getComposite(SMap, SI, " 153 -154 ");
	  //	  else Out += ModelSupport::getComposite(SMap, SI, " 114 134 ");
	  addOuterUnionSurf(Out);
	  
	}
      }
      
    // Flight lines 
    // x<0
      Out = ModelSupport::getComposite(SMap, SI, " 105 -106 -5 (103 -104  ");
      if (FlightLineType==1) Out += ModelSupport::getComposite(SMap, SI, " : (((-103 -143) : (104 -143))  ") + ReflectorSideAl + ")";
      if (TopPreType) {
      } else {
	Out += ModelSupport::getComposite(SMap, SI, SI-50, ": ((9M 39M -1M -2M (2 : 12)) : (4 14 -2 -12))"); // need this when W=H=50
      }
      Out += ")" + ModComplement + " " + sSurf; // ) due to 105 -106
      System.addCell(MonteCarlo::Qhull(cellIndex++, 0, 0.0, Out));

      //      HeadRule HR;
      HR.procString(Out);
      HR.makeComplement();

      Out = ModelSupport::getComposite(SMap, SI, " 115 -116 -5 ( 113 -114 ");
      if (FlightLineType==1)	Out += ModelSupport::getComposite(SMap, SI, " : (((-113 -153) : (114 -153))  ") + ReflectorSideBe + ")";
      if (TopPreType) {
      } else {
	Out += ModelSupport::getComposite(SMap, SI, SI-50, ": ((9M 39M -1M -2M (2:12)) : (4 14 -2 -12))");
      }
      Out += ")" + ModComplement;
      Out1 = Out;
      Out +=  " ( " + HR.display() + " ) ";
      Out += ModelSupport::getComposite(SMap, SI-50, " 48M 46M ");
      System.addCell(MonteCarlo::Qhull(cellIndex++, FlightLineWallMat, 0.0, Out+sSurf));

      addOuterUnionSurf(Out1);
      

      // x>0
      Out = ModelSupport::getComposite(SMap, SI, " 125 -126 5 (123 -124 ");
      if (FlightLineType==1) Out += ModelSupport::getComposite(SMap, SI, " : (((-123 144) : (124 144)) ") + ReflectorSideAl + ")";
      if (TopPreType) {
      } else {
	Out += ModelSupport::getComposite(SMap, SI, SI-50, ": ((19M 29M 7M 8M (22:32)) : (24 34 -22 -32))"); // need this when W=H=50
      }
      Out += ")" + ModComplement + " " + sSurf;
      System.addCell(MonteCarlo::Qhull(cellIndex++, 0, 0.0, Out));

      HR.procString(Out);
      HR.makeComplement();
 
      Out = ModelSupport::getComposite(SMap, SI, " 135 -136 5 (133 -134 ");
      if (FlightLineType==1) Out += ModelSupport::getComposite(SMap, SI, " : (((-133 154) : (134 154)) ") + ReflectorSideBe + ")";
      if (TopPreType) {
      } else
	Out += ModelSupport::getComposite(SMap, SI, SI-50, ": ((19M 29M 7M 8M (22 : 32)) : (24 34 -22 -32))"); // !!! removed 125 -126 - is it correct?
      Out +=  ")" + ModComplement;
      Out1 = Out;
      Out += " ( " + HR.display() + " ) ";
      Out += ModelSupport::getComposite(SMap, SI-50, " 48M 46M ");
      System.addCell(MonteCarlo::Qhull(cellIndex++, FlightLineWallMat, 0.0, Out+sSurf));

      addOuterUnionSurf(Out1);

    
    return; 
  }


  void ButterflyModerator::createAll(Simulation& System, const attachSystem::FixedComp& FC)
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

    ELog::RegMethod RegA("ButterflyModerator","createAll(Simulation &, const attachSystem::FixedComp&)");
    ELog::EM << "!!! This method should not be called !!!" << ELog::endErr;

    return;
  }



  void ButterflyModerator::createAll(Simulation& System, const attachSystem::FixedComp& FC, const attachSystem::FixedComp& ShutterBay, const long int sIndex)
  {
    /*!
      Extrenal build everything
      \param System :: Simulation
      \param FC :: FixedComponent for origin

      In our case the reflector has to be built relative to an origin and an axes set. 
      If you take a simple fixed object, then the axes is the axes set of this fixed object and the origin is the origin of this fixed object,
      which does not mean that the reflector and the object have the same origin, that's just the way you start and then you add the next bits.
    */

    ELog::RegMethod RegA("ButterflyModerator","createAll");
    // the order matters:

    populate(System.getDataBase()); // populate variables
    createUnitVector(FC); // take fixed component, then apply shift and angle rotation (transformation) for this object centre

    int choice = 0;
    if (choice == 0) { // butterfly
      createSurfaces();
      createObjects(System, ShutterBay, sIndex);
    } else if (choice == 1) { // rectangle
      createSurfacesREC();
      createObjectsREC(System, ShutterBay, sIndex);
    }

    createLinks();
    insertObjects(System);       

    return;
  }

  void
  ButterflyModerator::createLinks()
  /*!
    Creates a full attachment set
    Links/directions going outwards true.
  */
  {
    ELog::RegMethod RegA("ButterflyModerator", "createLinks");

    const size_t NL = nLayers-1;
    const int SI(modIndex+static_cast<int>(NL)*50);

    FixedComp::setConnect(0, Origin-Y*Length[NL]/2.0, -Y);
    FixedComp::setLinkSurf(0, -SMap.realSurf(SI+3));
    FixedComp::addBridgeSurf(0, -SMap.realSurf(SI+11));
    
    FixedComp::setConnect(1, Origin+Y*Length[NL]/2.0, Y);
    FixedComp::setLinkSurf(1, SMap.realSurf(SI+4));
    FixedComp::addBridgeSurf(1, SMap.realSurf(SI+11));
    
    // 2,3 - x

    FixedComp::setConnect(4, Origin-Z*Height[NL]/2.0, -Z);
    FixedComp::setLinkSurf(4, -SMap.realSurf(SI+5));
    FixedComp::addBridgeSurf(4, -SMap.realSurf(SI+5));

    FixedComp::setConnect(5, Origin+Z*Height[NL]/2.0, Z);
    FixedComp::setLinkSurf(5, SMap.realSurf(SI+6));
    FixedComp::addBridgeSurf(5, SMap.realSurf(SI+6));


    // links to the inner layer - for OnionStructure
    FixedComp::setConnect(6, Origin-Z*Height[0]/2.0, -Z);
    FixedComp::setLinkSurf(6, -SMap.realSurf(modIndex+5));
    //    FixedComp::addBridgeSurf(6, SMap.realSurf(modIndex+5));

    FixedComp::setConnect(7, Origin+Z*Height[0]/2.0, Z);
    FixedComp::setLinkSurf(7, SMap.realSurf(modIndex+6));
    //    FixedComp::addBridgeSurf(7, SMap.realSurf(modIndex+6));

    // TopPre moderator
    FixedComp::setConnect(8, Origin+Z*(Height[NL]/2.0+TopPreHeight), Z);
    FixedComp::setLinkSurf(8, SMap.realSurf(SI+50+55));

    return;
  }


  // virtual functions needed due to inheritance from LayerComp:
  Geometry::Vec3D
  ButterflyModerator::getSurfacePoint(const size_t layerIndex,
				const size_t sideIndex) const
  /*!
    Given a side and a layer calculate the link point
    \param sideIndex :: Side [0-1]
    \param layerIndex :: layer, 0 is inner moderator [0-4]
    \return Surface point
  */
  {
    ELog::RegMethod RegA("ButterflyModerator","getSurfacePoint");

    const int nSides = 9;
    if (sideIndex>nSides-1) 
      throw ColErr::IndexError<size_t>(sideIndex,nSides,"sideIndex ");
    if (layerIndex>nLayers) 
      throw ColErr::IndexError<size_t>(layerIndex,nLayers,"layer");

    switch(sideIndex) // kbat: see CylModerator.cxx
      {
      case 0:
	return Origin-Y*(Length[layerIndex]/2.0);
      case 1:
	return Origin+Y*(Length[layerIndex]/2.0);
      case 2:
	return Origin-X*(Width[layerIndex]/2.0);
      case 3:
	return Origin+X*(Width[layerIndex]/2.0);
      case 4:
	return Origin-Z*(Height[layerIndex]/2.0);
      case 5:
	return Origin+Z*(Height[layerIndex]/2.0);
      case 8:
	return Origin+Z*(Height[nLayers-1]/2.0+TopPreHeight);
      }
      throw ColErr::IndexError<size_t>(sideIndex,nSides,"sideIndex ");
  }

  std::string
  ButterflyModerator::getLayerString(const size_t layerIndex,
			       const size_t sideIndex) const
  /*!
    Given a side and a layer calculate the link surf
    \param sideIndex :: Side [0-5]
    \param layerIndex :: layer, 0 is inner moderator [0-4]
    \return Surface string
  */
  {
    ELog::RegMethod RegA("ButterflyModerator","getLayerString");

    if (layerIndex>nLayers) 
      throw ColErr::IndexError<size_t>(layerIndex,nLayers,"layer");

    const int SI(modIndex+static_cast<int>(layerIndex)*50);
    std::ostringstream cx;
    switch(sideIndex)
      {
      case 0:
	//	std::cout << "LayerString side 0 layer " << layerIndex << " "  << ModelSupport::getComposite(SMap,SI, " -3 ") << std::endl;
	return ModelSupport::getComposite(SMap,SI, " -3 ");
      case 1:
	return ModelSupport::getComposite(SMap,SI, " 4 ");
      case 2:
	throw ColErr::IndexError<size_t>(sideIndex,2,"sideIndex "); // -x
      case 3:
	throw ColErr::IndexError<size_t>(sideIndex,3,"sideIndex "); // +x
      case 4:
	return ModelSupport::getComposite(SMap, SI, " -5 ");
      case 5:
	return ModelSupport::getComposite(SMap, SI, " 6 ");
      }
    throw ColErr::IndexError<size_t>(sideIndex,6,"sideIndex ");
  }

  int
  ButterflyModerator::getLayerSurf(const size_t layerIndex,
			     const size_t sideIndex) const
  /*!
    Given a side and a layer calculate the link surf
    \param sideIndex :: Side [0-5]
    \param layerIndex :: layer, 0 is inner moderator [0-4]
    \return Surface string
  */
  {
    ELog::RegMethod RegA("ButterflyModerator","getLayerSurf");
    ELog::EM << "Not implemented " << ELog::endErr;

    if (layerIndex>nLayers) 
      throw ColErr::IndexError<size_t>(layerIndex,nLayers,"layer");
  
    const int SI(modIndex+static_cast<int>(layerIndex)*50);
    switch(sideIndex)
      {
      case 0:
	return -SMap.realSurf(SI+3);
      case 1:
	return  SMap.realSurf(SI+4);
      }
      throw ColErr::IndexError<size_t>(sideIndex,2,"sideIndex ");
      return 0;
  }

  int
  ButterflyModerator::getCommonSurf(const size_t sideIndex) const
  /*!
    Given a side calculate the boundary surface
    \param sideIndex :: Side [0-5]
    \return Common dividing surface [outward pointing]
  */
  {
    ELog::RegMethod RegA("ButterflyModerator","getCommonSurf");
    const int SI(modIndex+static_cast<int>(nLayers-1)*50);

    switch(sideIndex)
      {
      case 0:
	return -SMap.realSurf(SI+11);
      case 1:
	return SMap.realSurf(SI+11);
      case 2:
	return SMap.realSurf(SI+50+5);
      case 3:
	return SMap.realSurf(SI+50+5);
      case 4:
	return -SMap.realSurf(SI+12);
      case 5:
	return  SMap.realSurf(SI+12);
      case 8:
	return  SMap.realSurf(modIndex+50+55);
      }
      throw ColErr::IndexError<size_t>(sideIndex,6,"sideIndex ");
  }

  ButterflyModerator*
  ButterflyModerator::clone() const
  {
    /*!
      Clone copy constructor
      \return copy of this
    */
    return new ButterflyModerator(*this);
  }

  void ButterflyModerator::replaceSurface(std::string &s, int SIold, int SInew, const char *iold, const char *inew)
  {
    std::string sold = ModelSupport::getComposite(SMap, SIold, iold);
    std::string snew = ModelSupport::getComposite(SMap, SInew, inew);
    boost::replace_all(s, sold, snew);
  }

  std::string ButterflyModerator::getLinkComplement(const size_t sideIndex) const
  {
    return FixedComp::getLinkComplement(sideIndex);
  }
  
  std::string ButterflyModerator::getLinkString(const size_t sideIndex) const
  {
    return FixedComp::getLinkString(sideIndex);
  }


  //                                                 test of piping
  void ButterflyModerator::createSurfacesREC()  
  {
    // create surfaces for piping test
    int SI(modIndex);
    for (size_t i=0; i<nLayers; i++) {
      //      ModelSupport::buildPlane(SMap, SI+1, Origin-X*(Width[i]/2.0), X);
      ModelSupport::buildPlane(SMap, SI+1, Origin+Geometry::Vec3D(-Width[i]/2.0, 0, 0), X);
      ModelSupport::buildPlane(SMap, SI+2, Origin+X*(Width[i]/2.0), X);
      ModelSupport::buildPlane(SMap, SI+3, Origin-Y*(Length[i]/2.0), Y);
      ModelSupport::buildPlane(SMap, SI+4, Origin+Y*(Length[i]/2.0), Y);
      ModelSupport::buildPlane(SMap, SI+5, Origin-Z*(Height[i]/2.0), Z);
      ModelSupport::buildPlane(SMap, SI+6, Origin+Z*(Height[i]/2.0), Z);

      ModelSupport::buildPlane(SMap, SI+11, Origin, Y); // plane along the x-axis to separate y- from y+

      SI += 50;
    }
    ModelSupport::buildPlane(SMap, SI+5, Origin, X);
  }

  void ButterflyModerator::createObjectsREC(Simulation& System,
				    const attachSystem::FixedComp& ShutterBay, const long int sIndex)
  {
    std::string Out;
    int SI(modIndex);
    
    for (size_t i=0; i<nLayers; i++) {
      Out = ModelSupport::getComposite(SMap, SI, " 1 -2 3 -4 5 -6 ");
      if (i == nLayers-1) addOuterUnionSurf(Out);
      
      if (i) {
        Out += ModelSupport::getComposite(SMap, SI-50, " (-1:2:-3:4:-5:6) ");
      }

      System.addCell(MonteCarlo::Qhull(cellIndex++, mat[i], temp[i], Out));
      SI += 50;
    }

  }


  void ButterflyModerator::createCYLSurfaces()
  {
    ELog::RegMethod RegA("ButterflyModerator","createCYLSurfaces");
    
    // Divide plane
    ModelSupport::buildPlane(SMap,modIndex+1,Origin,X);  
    ModelSupport::buildPlane(SMap,modIndex+2,Origin,Y);  
    
    int SI(modIndex);
    for(size_t i=0;i<nLayers;i++)
      {
	ModelSupport::buildCylinder(SMap,SI+7,Origin,Z,Width[i]);  
	ModelSupport::buildPlane(SMap,SI+5,Origin-Z*Height[i]/2.0,Z);  
	ModelSupport::buildPlane(SMap,SI+6,Origin+Z*Height[i]/2.0,Z);  
	SI+=10;
      }
    
    return; 
  }

  void ButterflyModerator::createCYLObjects(Simulation& System)
  {
    ELog::RegMethod RegA("ButterflyModerator","createCYLObjects");

    std::string Out;
    int SI(modIndex);
    for(size_t i=0;i<nLayers;i++) {
      Out=ModelSupport::getComposite(SMap,SI," -7 5 -6 ");
    
      if ((i+1)==nLayers) addOuterSurf(Out);
      if (i)
	Out+=ModelSupport::getComposite(SMap,SI-10," (7:-5:6) ");
      System.addCell(MonteCarlo::Qhull(cellIndex++,mat[i],temp[i],Out));
      SI+=10;
    }
    return; 
  }




}  // NAMESPACE instrumentSystem
