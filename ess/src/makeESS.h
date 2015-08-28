/********************************************************************* 
  CombLayer : MNCPX Input builder
 
 * File:   essBuildInc/makeESS.h
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
#ifndef essSystem_makeESS_h
#define essSystem_makeESS_h

#include <boost/shared_array.hpp>

#define NF5 2

namespace constructSystem
{
  class ModBase;
  class SupplyPipe;
}

namespace essSystem
{
  class WheelBase;
  class SegWheel;
  class Wheel;
  class BeRef;
  class CylinderCell;
  class essMod;
  class ConicModerator;
  class CylModerator;
  class WaterPipe;
  class CylPreMod;
  class BulkModule;
  class ShutterBay;
  class GuideBay;
  class F5Collimator;
  class FlightLine;
  class ProtonVoid;
  class ProtonBeamWindow;
  class RecModerator;
  class ButterflyModerator;
  class ThermalModerator;
  class FlowGuide;
  class Grooving;

  /*!
    \class makeESS
    \version 1.0
    \author S. Ansell
    \date January 2013
    \brief WheeModerator for ESS 
  */
  
class makeESS
{
 private:
  
  std::shared_ptr<WheelBase> Target;   ///< target object
  std::shared_ptr<BeRef> TopReflector;   ///< top reflector object
  std::shared_ptr<BeRef> LowReflector;   ///< low reflector object

  std::shared_ptr<RecModerator> TubeMod;  ///< TubeModerator
  std::shared_ptr<Grooving> TubeModGrooving; ///< Tube moderator groove
  std::shared_ptr<ThermalModerator> ThermalRecMod;  ///rec thermal moderator

  std::shared_ptr<ButterflyModerator> LowButterfly;  ///< low Butterfly Modeartor
  std::shared_ptr<ButterflyModerator> TopButterfly;  ///< top Butterfly Modeartor

  std::shared_ptr<FlowGuide> OnionModPipe;   ///< onion-like cooling pipes for the flat moderator
  std::shared_ptr<FlowGuide> OnionPrePipe;   ///< onion-like cooling pipes for Premoderator
  std::shared_ptr<FlowGuide> LowModOnion;   ///< onion-like cooling pipes for the low moderator
  std::shared_ptr<FlowGuide> LowModLFlowGuide, LowModRFlowGuide, TopModLFlowGuide, TopModRFlowGuide;   ///< flow guide for top moderator - hydrogen part, left and right
  std::shared_ptr<essSystem::ProtonVoid> PBeam;  ///< Proton Void
  std::shared_ptr<essSystem::ProtonBeamWindow> PBW;  ///< Proton beam window


  // ASSEMBLY 1:
  std::shared_ptr<essSystem::CylModerator> LowMod;  ///< Lower Mod
  //  std::shared_ptr<constructSystem::ModBase> LowModB;  ///< Lower Mod [if needed]
  std::shared_ptr<essSystem::FlightLine> LowAFL;  ///< Lower Mode FL
  std::shared_ptr<essSystem::FlightLine> LowBFL;  ///< Lower Mode FL
  std::shared_ptr<CylPreMod> LowPre;  ///< Upper Mod (Pre)
  std::shared_ptr<constructSystem::SupplyPipe> LowSupplyPipe;  ///< Lower supply 
  std::shared_ptr<constructSystem::SupplyPipe> LowReturnPipe;  ///< Lower return
  std::shared_ptr<constructSystem::SupplyPipe> TopSupplyPipe;  ///< Upper supply 
  std::shared_ptr<constructSystem::SupplyPipe> TopReturnPipe;  ///< Upper return 

  std::shared_ptr<WaterPipe> TopPreAPipe;  ///< Upper pipe A of water premoderator
  std::shared_ptr<WaterPipe> TopPreBPipe;  ///< Upper pipe B of water premoderator

  std::shared_ptr<essSystem::CylModerator> TopMod;   ///< Upper Mod
  std::shared_ptr<essSystem::FlightLine> TopAFL;  ///< Upper Mode FL
  std::shared_ptr<essSystem::FlightLine> TopBFL;  ///< Upper Mode FL
  std::shared_ptr<CylPreMod> TopPre;  ///< Upper Mod (Pre)

  std::shared_ptr<BulkModule> Bulk;      ///< Main bulk module
  std::shared_ptr<essSystem::FlightLine> BulkLowAFL;  ///< Lower Mode FL (kbat: maybe not used?)

  // cylinder of water between the target wheel and the bottom thermal moderator. It's of the BeRef type since this class is a good choice for a cylinder.
  std::shared_ptr<BeRef> LowWaterDisc;   ///< low water disc
  std::shared_ptr<BeRef> TopWaterDisc;   ///< top water disc
  std::shared_ptr<CylinderCell> ThermalCylMod; ///< thermal moderator in the configuration of cylindrical thermal moderator
  

  /// Shutterbay object
  std::shared_ptr<ShutterBay> ShutterBayObj;  
  std::vector<std::shared_ptr<GuideBay> > GBArray;  
  
  // F5 collimators
  //  int nF5col; // number of F5 collimators
  std::vector<std::shared_ptr<F5Collimator>> F5array;



  void topFlightLines(Simulation&, const std::string lowModType);
  void lowFlightLines(Simulation&);
  void createGuides(Simulation&, const std::string lowModType);
  void buildLowMod(Simulation&);
  void buildF5Collimator(Simulation&, std::shared_ptr<F5Collimator> F5);
  //  void buildConicMod(Simulation&);
  void makeTarget(Simulation&,const mainSystem::inputParam&);

  void dumpMaterialMesh(Simulation& SimPtr, const Geometry::Vec3D &startPt, const Geometry::Vec3D &endPt, const size_t nX, const size_t nY, const size_t nZ, const char *fname) const;
  std::vector<int> getMaterials(Simulation& SimPtr, const Geometry::Vec3D &center, double *stepXYZ, size_t N) const;
  std::string getMaterialString(std::vector<int>) const;

 public:
  
  makeESS();
  makeESS(const makeESS&);
  makeESS& operator=(const makeESS&);
  ~makeESS();
  
  void build(Simulation*,const mainSystem::inputParam&);

};

}

#endif
