#ifndef essSystem_RecModerator_h
#define essSystem_RecModerator_h

class Simulation;

namespace essSystem
{

  /*!
    \class RecModerator
    \author S. Ansell
    \version 1.0
    \date February 2013
    \brief Reflector object 
  */

  class RecModerator : public constructSystem::ModBase
      {
      private:
	
	std::vector<double> Width; ///< along x
	std::vector<double> Depth; ///< along y
	std::vector<double> Height; ///< along z
	std::vector<int> mat; ///< materials
	std::vector<double> temp; ///< temperatures

	double windowThick; ///< thickness of Al window where neutrons are extracted. Must be less than Width[nLayers-1]

	double preWidth;                /// water premoderator thickness along x
	double preDepth;                /// along y
	double preHeight1;               /// along z - lower water layer
	double preHeight2;               /// along z - upper water layer
	double preWingThick;               /// thickness of permoderator wings

	int preWallMat;                    ///< wall Material
	double preWallThick;               ///< wall thickness
	int preMat;                     ///< premoderator material (H2O)
	double preTemp;                 ///< premoderator temperature

	double tubeDepth;    ///< extraction tube depth (along y)
	double tubeHeight;    ///< extraction tube height (along z)

	int targetGapMat; ///< material of the gap between the target wheel and premoderator
	int fastRefMat;  ///< material of fast reflector

	double WaterCylinderRad; //< radius of the water cylinder above the rectangular moderator
	double WaterCylinderMat; //< material of the water cylinder above the rectangular moderator

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

	RecModerator(const std::string&);
	RecModerator(const RecModerator&);
	RecModerator& operator=(const RecModerator&);
	virtual ~RecModerator();

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
	virtual RecModerator* clone() const;

      };

}

#endif
 
