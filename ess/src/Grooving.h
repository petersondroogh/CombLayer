#ifndef essSystem_Grooving_h
#define essSystem_Grooving_h

class Simulation;

namespace essSystem
{

  /*!
    \class Grooving
    \author Konstantin Batkov
    \version 1.0
    \date July 2014
    \brief Grooving object
  */

  class Grooving : public attachSystem::ContainedComp, public attachSystem::FixedComp
    {
    private:

      const int grIndex;             ///< Index of surface offset
      int cellIndex;                  ///< Cell index. Every object has a space of about 10000 cells unless you request more.

      double xStep;                   ///< X step
      double yStep;                   ///< Y step
      double zStep;                   ///< Z step
      double xyAngle;                 ///< XY Angle
      double zAngle;                  ///< Z Angle

      double length;
      double height;
      double width;

      double wallThick;               ///< Wall thickness

      int nGrooves;                   ///< number of grooves (along y)
      double dy;                      ///< distance between outer surfaces of neighbouring grooves
      int grMat;                      ///< groove material - usually void
      int wallMat;                    ///< wall Material

      // Functions:

      void populate(const FuncDataBase&);
      void createUnitVector(const attachSystem::FixedComp&);

      void createSurfaces();
      void createObjects(Simulation&, const std::string& s1, const std::string& s2);
      void createLinks();

    public:

      Grooving(const std::string&);
      Grooving(const Grooving&);
      Grooving& operator=(const Grooving&);
      virtual ~Grooving();

      virtual void addToInsertChain(attachSystem::ContainedComp&) const;
      void createAll(Simulation&,const attachSystem::FixedComp&, const attachSystem::FixedComp&, const long int l1, const long int l2);

      inline int getGrIndex() { return grIndex; }
      inline int getCellIndex() { return cellIndex; }
  
    };

}

#endif
 
