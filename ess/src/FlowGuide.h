#ifndef essSystem_FlowGuide_h
#define essSystem_FlowGuide_h

class Simulation;

namespace essSystem
{

/*!
  \class FlowGuide
  \date May 2015
  \brief Flow guide structures
*/

class FlowGuide : public attachSystem::ContainedComp,
    public attachSystem::FixedComp
{
 private:

  const int onionIndex;             ///< Index of surface offset
  int cellIndex;                  ///< Cell index. Every object has a space of about 10000 cells unless you request more.

  double xStep;                   ///< X step
  double yStep;                   ///< Y step
  double zStep;                   ///< Z step
  double xyAngle;                 ///< XY Angle
  double zAngle;                  ///< Z Angle

  double wallThick;               ///< Wall thickness
  int wallMat;                        ///< material
  double wallTemp;  ///< wall temperature

  std::string Type; /// Flow guide type (Onion, Star)

  // variables for type 'Onion'
  size_t nRings;
  std::vector<double> radius;                  ///< Radius of the rings
  std::vector<double> gateWidth;                  ///< full width of spacing in the corresponding ring
  std::vector<double> gateLength;                  ///< length of 'the door' in the corresponding ring

  std::string BottomSurface; ///< bottom surface number
  std::string UpperSurface;  ///< upper surface number

  // variables for type 'Star'
  /*
 base  side
   |   /        forearm
   |  /         /
   ------- arm /
   |  \
   |   \
 base   side
  */
  double BaseYShift;
  double BaseLength;
  double BaseArmDist;
  double ArmYShift;
  double ArmLength; // length of one half
  double ForeArmLength;
  double ForeArmAngle;
  double SideXShift;
  double SideYShift;
  double SideLength;
  double SideAngle;

  // Functions:

  void populate(const FuncDataBase&);
  void createUnitVector(const attachSystem::FixedComp&);

  void createSurfaces();
  void createObjects(Simulation&);
  void createLinks();

 public:

  FlowGuide(const std::string&);
  FlowGuide(const FlowGuide&);
  FlowGuide& operator=(const FlowGuide&);
  virtual ~FlowGuide();

  void setBottomSurface(const attachSystem::FixedComp& FC, const long int link);
  void setUpperSurface(const attachSystem::FixedComp& FC, const long int link);

  int getMainCell() const { return onionIndex+1; }
  virtual void addToInsertChain(attachSystem::ContainedComp&) const; // has to be there [2:1040]
  void createAll(Simulation&,const attachSystem::FixedComp&); // [2:1070]
  
};

}

#endif
 
