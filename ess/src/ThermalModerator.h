#ifndef essSystem_ThermalModerator_h
#define essSystem_ThermalModerator_h

class Simulation;

namespace essSystem
{

  /*!
    \class ThermalModerator
    \author S. Ansell
    \version 1.0
    \date February 2013
    \brief Reflector object 
  */

  class ThermalModerator : public constructSystem::ModBase
      {
      private:
	
	std::vector<double> Width; ///< along x
	std::vector<double> Depth; ///< along y
	std::vector<double> Height; ///< along z
	std::vector<int> mat; ///< materials
	std::vector<double> temp; ///< temperatures

	double WaterCylinderRad; //< radius of the water cylinder above the rectangular moderator
	double WaterCylinderMat; //< material of the water cylinder above the rectangular moderator

	double AFlightAngleZBase, AFlightAngleZTop; //< x<0
	double BFlightAngleZBase, BFlightAngleZTop; //< x>0
	double AFlightAngleXY1, AFlightAngleXY2;
	double BFlightAngleXY1, BFlightAngleXY2;
	double FlightLineWallThick; ///< thickness of flight line wall
	int    FlightLineWallMat; ///< material of flight line wall


	// Functions:

	void populate(const FuncDataBase&);
	void createUnitVector(const attachSystem::FixedComp&);

	void createSurfaces();
	void createObjects(Simulation&,
			   const attachSystem::FixedComp& Target, const long int tIndex,
			   const attachSystem::FixedComp& ShutterBay, const long int sIndex,
			   const attachSystem::FixedComp& Reflector, const long int rSide, const long int rTop);
	void createLinks();

      public:

	ThermalModerator(const std::string&);
	ThermalModerator(const ThermalModerator&);
	ThermalModerator& operator=(const ThermalModerator&);
	virtual ~ThermalModerator();

	int getMainCell() const { return modIndex+1; }
	virtual void addToInsertChain(attachSystem::ContainedComp&) const; // has to be there [2:1040]
	void createAll(Simulation&,const attachSystem::FixedComp&); // [2:1070]
	void createAll(Simulation&,const attachSystem::FixedComp& FC, 
		       const attachSystem::FixedComp& Target, const long int tIndex,
		       const attachSystem::FixedComp& ShutterBay, const long int sIndex,
		       const attachSystem::FixedComp& Reflector, const long int rSide, const long int rTop); // [2:1070]
  

	virtual Geometry::Vec3D getSurfacePoint(const size_t,const size_t) const;
	virtual int getLayerSurf(const size_t,const size_t) const;
	virtual std::string getLayerString(const size_t,const size_t) const;
	virtual int getCommonSurf(const size_t) const;
	virtual ThermalModerator* clone() const;

      };

}

#endif
 
