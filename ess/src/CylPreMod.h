#ifndef essSystem_CylPreMod_h
#define essSystem_CylPreMod_h

class Simulation;

namespace essSystem
{
  class BlockAddition;
  /*!
    \class CylPreMod
    \author S. Ansell
    \version 1.0
    \date October 2012
    \brief Specialized for a cylinder moderator
  */

  class CylPreMod : public attachSystem::ContainedGroup,
    public attachSystem::LayerComp,
    public attachSystem::FixedComp
      {
      private:
  
	const int modIndex;             ///< Index of surface offset
	int cellIndex;                  ///< Cell index

	double xStep;                   ///< X step
	double yStep;                   ///< Y step
	double zStep;                   ///< Z step
	double xyAngle;                 ///< XY Angle
	double zAngle;                  ///< Z Angle
	
	/// Extension object
	std::shared_ptr<BlockAddition> ExtAObj; 
	std::shared_ptr<BlockAddition> ExtBObj; 
	int blockActiveA;
	int blockActiveB;
	size_t aSide;  
	size_t bSide;

	double innerRadius;             ///< Radius from inner cell
	double innerHeight;             ///< height from inner cell
	double innerDepth;              ///< Depth from inner cell

	std::vector<double> radius;         ///< cylinder radii
	std::vector<double> height;         ///< Full heights
	std::vector<double> depth;          ///< full depths
	std::vector<int> mat;               ///< Materials
	std::vector<double> temp;           ///< Temperatures

	std::vector<Geometry::Vec3D> viewY; ///< Direction from centre
	std::vector<double> viewAngle;      ///< Angle from centre
	std::vector<double> viewOpenAngle;  ///< Angle opening
	std::vector<double> viewHeight;     ///< Height from centre
	std::vector<double> viewWidth;      ///< Width at intercept

	// Now calculated
	std::vector<Geometry::Vec3D> FLpts;   ///< Flight line corner 
	std::vector<Geometry::Vec3D> FLunit;  ///< Flight line direction  [-x,x,-z,z]
	std::vector<int> layerCells;          ///< Layered cells
	// Functions:  
	void checkItems(const attachSystem::FixedComp&);
  
	void populate(const FuncDataBase&);
	void createUnitVector(const attachSystem::FixedComp&);

	void createSurfaces();
	void createObjects(Simulation&,const attachSystem::ContainedComp*, const attachSystem::FixedComp&, const long int rIndex);
	void createLinks();
  
	void updateLayers(Simulation&,const char,
			  const size_t,const size_t) const;

	Geometry::Vec3D calcViewIntercept(const size_t,const size_t) const;

      public:

	CylPreMod(const std::string&);
	CylPreMod(const CylPreMod&);
	CylPreMod& operator=(const CylPreMod&);
	virtual ~CylPreMod();

	virtual void addToInsertChain(attachSystem::ContainedComp&) const;
	void createAll(Simulation&,const attachSystem::FixedComp&, const attachSystem::FixedComp&, const long int rIndex);

	const std::shared_ptr<BlockAddition>& getBox(const char) const;
	std::string getBoxCut(const char) const;
	virtual Geometry::Vec3D getSurfacePoint(const size_t,const size_t) const;
	virtual int getLayerSurf(const size_t,const size_t) const;
	virtual std::string getLayerString(const size_t,const size_t) const;


      };

}

#endif
 
