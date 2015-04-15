#ifndef essSystem_FlightLine_h
#define essSystem_FlightLine_h

class Simulation;

namespace essSystem
{

  /*!
    \class FlightLine
    \version 1.0
    \author S. Ansell
    \date April 2011
    \brief FlightLine [insert object]
  */

  class FlightLine : public attachSystem::ContainedGroup,
    public attachSystem::FixedComp
      {
      private:
  
	const int flightIndex;        ///< Index of surface offset
	int cellIndex;                ///< Cell index
  
	double xStep;                 ///< Offset on X to Target
	double zStep;                 ///< Offset on Z top Target

	double masterXY;              ///< Master rotation of general axis(XY)
	double masterZ;               ///< Master rotation of general axis(Z)
	double anglesXY[2];           ///< Rotation in the XY plane 
	double anglesZ[2];            ///< Rotation in the Z plane
  
	double height;                ///< Height of flight line
	double width;                 ///< Width of flight line

	int plateIndex;               ///< Index of the side [+1] with sign 
	size_t nLayer;                ///< Number of layers
	std::vector<double> lThick;   ///< Liner Thickness 
	std::vector<double> lHeight;   ///< Liner Height
	std::vector<int> lMat;        ///< Layer Material

	bool capActive;
	std::vector<int> capLayer;    ///< End cap layers
	std::vector<HeadRule> capRule;   ///< Rule for each cap

	std::string attachRule;       ///< Attached rule
  
	void populate(const Simulation&);
	void createUnitVector(const attachSystem::FixedComp&,const size_t);
	void createUnitVector(const attachSystem::FixedComp&,const size_t, const size_t);
	void createUnitVector(const Geometry::Vec3D&,const Geometry::Vec3D&, const Geometry::Vec3D&);
	void createRotatedUnitVector(const attachSystem::FixedComp&,const size_t, const size_t);

	void createSurfaces();
	void createCapSurfaces(const attachSystem::FixedComp&,const size_t);
	void createObjects(Simulation&,const attachSystem::FixedComp&, const size_t);
	void createObjects(Simulation&,const attachSystem::FixedComp&, const int,const size_t, const attachSystem::ContainedComp&);

	void removeObjects(Simulation&);
	std::string getRotatedDivider(const attachSystem::FixedComp&, const size_t);

      public:

	FlightLine(const std::string&);
	FlightLine(const FlightLine&);
	FlightLine& operator=(const FlightLine&);
	~FlightLine();

	inline int getFlightIndex() const {return flightIndex;}
	inline int getCellIndex() const {return cellIndex;}

	void getInnerVec(std::vector<int>&) const;

	void createAll(Simulation&,const size_t, const attachSystem::FixedComp&);
	void createAll(Simulation&,const size_t,const size_t, const attachSystem::FixedComp&);
	void createAll(Simulation&,const attachSystem::FixedComp&, const attachSystem::ContainedComp&);

	void reBoundary(Simulation&,const size_t, const attachSystem::FixedComp&);

	void processIntersectMajor(Simulation&, const attachSystem::ContainedGroup&, const std::string&,const std::string&);
	void processIntersectMinor(Simulation&, const attachSystem::ContainedGroup&, const std::string&,const std::string&);
      };

}

#endif
 
