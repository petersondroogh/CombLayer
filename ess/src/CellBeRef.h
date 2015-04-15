#ifndef essSystem_CellBeRef_h
#define essSystem_CellBeRef_h

class Simulation;

namespace essSystem
{
  
  /*!
    \class CellBeRef
    \version 1.0
    \date May 2014
    \brief Reflector object with cell cooling
  */
  
  class CellBeRef : public attachSystem::ContainedComp, public attachSystem::FixedComp
    {
    private:
      
      const int refIndex;             ///< Index of surface offset
      int cellIndex;                  ///< Cell index. Every object has a space of about 10000 cells unless you request more.
      
      double xStep;                   ///< X step
      double yStep;                   ///< Y step
      double zStep;                   ///< Z step
      double xyAngle;                 ///< XY Angle
      double zAngle;                  ///< Z Angle

      double radius;                  ///< Radius
      double height;                  ///< Height
      double wallThick;               ///< Wall thickness

      int refMat;                     ///< reflector material - lower half
      //  int refMatTop;                     ///< reflector material - upper half
      int wallMat;                    ///< wall Material

      // Functions:

      void populate(const FuncDataBase&);
      void createUnitVector(const attachSystem::FixedComp&);

      void createSurfaces();
      void createObjects(Simulation&, const std::string&);
      void createLinks();

    public:

      CellBeRef(const std::string&);
      CellBeRef(const CellBeRef&);
      CellBeRef& operator=(const CellBeRef&);
      virtual ~CellBeRef();

      int getMainCell() const { return refIndex+1; } // Be
      int getWallCell() const { return refIndex+2; } // wall
      virtual void addToInsertChain(attachSystem::ContainedComp&) const; // has to be there [2:1040]
      void createAll(Simulation&,const attachSystem::FixedComp&, const attachSystem::FixedComp&, const long int); // [2:1070]
  
    };

}

#endif
 
