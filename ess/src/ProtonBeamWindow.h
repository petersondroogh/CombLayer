#ifndef essSystem_ProtonBeamWindow_h
#define essSystem_ProtonBeamWindow_h

class Simulation;

namespace essSystem
{
  class ProtonBeamWindow : public attachSystem::ContainedComp, public attachSystem::FixedComp {
      private:
	const int pbwIndex;             ///< Index of surface offset
	int cellIndex;                  ///< Cell index

	bool active;                    /// true if active
	double xStep;                   ///< X step (point detector x coordinate)
	double yStep;                   ///< Y step (point detector y coordinate)
	double zStep;                   ///< Z step (point detector z coordinate)
	double xyAngle;                 ///< XY Angle
	double zAngle;                  ///< Z Angle

	int linacMat;                   ///< material from the accelerator side
	int monolithMat;                ///< material from the monolith side

	double pipeRadius;               /// pipe radius
	double pipeWallThick;            /// pipe wall thickness
	int pipeWallMat;              /// pipe wall material
	int pipeCoolingMat;           /// pipe cooling material
	int nPipes;                    /// number of pipes

	double frameHeight;
	double frameWidth;
	double frameThick;
	int    frameMat;

	double panpipeHeight;
	double panpipeWidth;


	// Functions:

	void populate(const FuncDataBase&);
	void createUnitVector(const attachSystem::FixedComp&);

	void createSurfaces(const attachSystem::FixedComp&);
	void createObjects(Simulation&, const std::string&);
	void createLinks();


      public:
	ProtonBeamWindow(const std::string&);
	ProtonBeamWindow(const ProtonBeamWindow&);
	ProtonBeamWindow& operator=(const ProtonBeamWindow&);
	virtual ~ProtonBeamWindow();

	int getMainCell() const { return pbwIndex+1; }
	virtual void addToInsertChain(attachSystem::ContainedComp&) const; 
	void createAll(Simulation&,const attachSystem::FixedComp&, const long int);
  };
}

#endif
