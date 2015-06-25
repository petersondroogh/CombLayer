#ifndef essSystem_ButterflyModerator_h
#define essSystem_ButterflyModerator_h

class Simulation;

namespace essSystem
{


  /*!
    \class ButterflyModerator
    \author Konstantin Batkov
    \version 1.0
    \date December 2014
    \brief Butterfly moderator obect - implementation of the Troels' idea
  */

  class ButterflyModerator : public constructSystem::ModBase
      {
      private:
	
	std::vector<double> Width; ///< along x
	std::vector<double> Length; ///< along y
	std::vector<double> Height; ///< along z
	std::vector<int> mat; ///< materials
	std::vector<double> temp; ///< temperatures

	double WingLength; //< wing length (along y) of the Hydrogen part
	double WingCurvatureSideRadius; // curvature radius of the side wing edges of the Hydrogen part
	double WingCurvatureCentralRadius; // curvature radius of the central wing edge of the Hydrogen part

	double AFlightAngleZBase, AFlightAngleZTop; //< x<0
	double BFlightAngleZBase, BFlightAngleZTop; //< x>0
	double AFlightAngleXY1, AFlightAngleXY2;
	double BFlightAngleXY1, BFlightAngleXY2;
	double FlightLineWallThick; ///< thickness of flight line wall
	int    FlightLineWallMat; ///< material of flight line wall
	double FlightXOffset; ///< distance along the x-axis between the end of the wing and the start of the flight line plane
	double FlightTiltOffset; ///< coordinate along the x-axis where the vertical tilting starts
	int    FlightLineType; ///< 0=normal flight line; 1=planes inside BeRef are parallel to the x-axis
	int    FlightLineWrapTopPreType; //< Describes how flight line wraps TopPremoderator. 0=TopPre inside void of flight line; 1=TopPre inside inner reflector;
	
	// thermal moderator (side premoderator)
	int    PreMat;
	double PreTemp;
	double PreLength;
	double PreThick1; // water thickness outside para hydrogen
	double PreThick2; // water thickness inside (under/above) para hydrogen
	double PreBaseThick; // outside thickness at y=0
	double PreHeight;
	double PreWallThick;
	int    PreWallMat;
	double    PreWallTemp;
	double PreTopHeight; // thickness of water layer above cold moderator
	double PreBottomHeight; // thickness of water layer below cold moderator

	bool WingsSeparated; // true if the butterfly wings are separated, otherwise false
	int PipeType; // type of supply/return pipes (1=horizontal, 2=vertical. It's not enough just to change PipeType - the PPt* vectors must be adjusted accordingly)

	int TopPreType; // 0=no TopPre, 1=covers Butterfly, 2=covers BeRef at butterfly width; 3=cylinder of radius TopPreRadius
	double TopPreRadius; // makes sence only if TopPreType==3
	double TopPreHeight; // top premoderator (water layer)
	double TopPreWidth; // along x, makes sense only if TopPreType==2. It does not make sense to set it > BeRefRad, because then it becomes a cylinder
	int TopPreCoolingChannels; // 0=no cooling channels, 1=yes
	double TopPreTopVoidHeight; ///< height of void layer above TopPre. Implemented for TopPreType==1
	double TopPreTopAlHeight;   ///< height of the second Al layer above TopPre (and other Butterfly parts). Implemented for TopPreType==1

	int MaterialInBeRef; // material inside BeRef at the level of the moderator (usually water)

	std::string ReflectorSideAl; // BeRef side surface + Al container
	std::string ReflectorSideBe; // BeRef side surface

	// Functions:

	void populate(const FuncDataBase&);
	void createUnitVector(const attachSystem::FixedComp&);

	void createSurfaces();
	void createObjects(Simulation&, const attachSystem::FixedComp& ShutterBay, const long int sIndex);

	// methods to test a simple rectangular layered geometry
	void createSurfacesREC();
	void createObjectsREC(Simulation&, const attachSystem::FixedComp& ShutterBay, const long int sIndex);

	void createCYLSurfaces();
	void createCYLObjects(Simulation&);


	void createLinks();

	void replaceSurface(std::string &s, int SIold, int SInew, const char *iold, const char *inew);

      public:

	ButterflyModerator(const std::string&);
	ButterflyModerator(const ButterflyModerator&);
	ButterflyModerator& operator=(const ButterflyModerator&);
	virtual ~ButterflyModerator();

	int getMainCell() const { return modIndex+1; }
	virtual void addToInsertChain(attachSystem::ContainedComp&) const; // has to be there [2:1040]
	void createAll(Simulation&,const attachSystem::FixedComp&); // [2:1070]
	void createAll(Simulation&,const attachSystem::FixedComp& FC, 
		       const attachSystem::FixedComp& ShutterBay, const long int sIndex);

	void setReflectorSurfaces(const attachSystem::FixedComp& Reflector, const long int rSideAl, const long int rSideBe);
  

	virtual Geometry::Vec3D getSurfacePoint(const size_t,const size_t) const;
	virtual int getLayerSurf(const size_t,const size_t) const;
	virtual std::string getLayerString(const size_t,const size_t) const;
	virtual int getCommonSurf(const size_t) const;
	virtual ButterflyModerator* clone() const;

	virtual std::string getLinkComplement(const size_t) const;
	virtual std::string getLinkString(const size_t) const;

	int getPipeType() const { return PipeType; }
	inline bool IsTopPreCoolingChannels() const { return TopPreCoolingChannels; }
	inline int getTopPreType() const { return TopPreType; }
      };

}

#endif
 
