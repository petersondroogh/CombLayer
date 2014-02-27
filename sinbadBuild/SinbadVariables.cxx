#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>
#include <list>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <algorithm>
#include <iterator>
#include <boost/array.hpp>
#include <boost/shared_ptr.hpp>

#include "Exception.h"
#include "FileReport.h"
#include "NameStack.h"
#include "RegMethod.h"
#include "GTKreport.h"
#include "OutputLog.h"
#include "support.h"
#include "MatrixBase.h"
#include "Matrix.h"
#include "Vec3D.h"
#include "Code.h"
#include "varList.h"
#include "FuncDataBase.h" 
#include "variableSetup.h"

namespace setVariable
{

void
SinbadVariables(FuncDataBase& Control)
  /*!
    Function to set the control variables and constants
    -- This version is for Sinbad ()
    \param Control :: Function data base to add constants too
  */
{
  // Detectors

  Control.addVariable("49DetectorPositionN",4);  

  Control.addVariable("49Detector0Active",1);   
  Control.addVariable("49Detector1Active",0);   
  Control.addVariable("49Detector2Active",0);   
  Control.addVariable("49Detector3Active",0);   

  Control.addVariable("49DetectorXStep",2.5);  
  Control.addVariable("49DetectorYStep",0.0);  
  Control.addVariable("49DetectorZStep",-6.6);  
  Control.addVariable("49DetectorLength","Void");  
  Control.addVariable("49DetectorMat","Rh");  

  // Rh
  Control.addVariable("49DetectorRadius",1.27);   
  Control.addVariable("49DetectorLength",0.1); 

  Control.addVariable("49DetectorYStep1",5.1);  
  Control.addVariable("49DetectorYStep2",10.18);  
  Control.addVariable("49DetectorYStep3",15.96);  
  Control.addVariable("49DetectorYStep4",21.74);  
  Control.addVariable("49DetectorYStep5",25.44);  
  Control.addVariable("49DetectorYStep6",27.44);  
  Control.addVariable("49DetectorYStep7",29.44);  

  Control.addVariable("49DetectorYStep8",29.94);  
  Control.addVariable("49DetectorYStep9",32.94);  
  Control.addVariable("49DetectorYStep10",34.44);  

  //Nestor side

  Control.addVariable("NestorXStep",91.45);  
  Control.addVariable("NestorYStep",-1.45);  
  Control.addVariable("NestorZStep",95.5);  

  Control.addVariable("NestorWidth",182.9);   
  Control.addVariable("NestorHeight",191.0);  

  Control.addVariable("NestorNSlab",6);
  Control.addVariable("NestorThick0",0.70); 
  Control.addVariable("NestorMat0","Void"); 
  Control.addVariable("NestorThick1",15.0); 
  Control.addVariable("NestorMat1","Graphite"); 
  Control.addVariable("NestorThick2",0.60); 
  Control.addVariable("NestorMat2","Void");  
  Control.addVariable("NestorThick3",5.08); 
  Control.addVariable("NestorMat3","sbadLead");  
  Control.addVariable("NestorThick4",0.70); 
  Control.addVariable("NestorMat4","Void");  
  Control.addVariable("NestorThick5",3.18); 
  Control.addVariable("NestorMat5","sbadMildSteel");  
  Control.addVariable("NestorAlWindowRadius",56.06); 








  // define centre
  Control.addVariable("49ShieldOffSetX",91.45);  
  Control.addVariable("49ShieldOffSetY",1.45);  
  Control.addVariable("49ShieldOffSetZ",95.5);  

  Control.addVariable("49ShieldLengthSlab",182.9);   
  Control.addVariable("49ShieldHeightSlab",191.);  
  Control.addVariable("49ShieldTemperatureSlab",300.0); 

  Control.addVariable("49ShieldFrontSlabN",26);
  
 // this void cell moved to accomodate detector after boral
  // Control.addVariable("49ShieldWidthSlab0",4.60); 

  // this for S 
  // Control.addVariable("49ShieldWidthSlab0",4.04);
  // this for Rh
  Control.addVariable("49ShieldWidthSlab0",4.5985); 
  // this for base Mn
    // Control.addVariable("49ShieldWidthSlab0",4.585); 
    // Mn+cd
    // Control.addVariable("49ShieldWidthSlab0",4.458); 
 
  Control.addVariable("49ShieldMaterialSlab0","Void"); 

  Control.addVariable("49ShieldWidthSlab1",0.50); 
  Control.addVariable("49ShieldMaterialSlab1","Boral5"); 
  // Rh
    Control.addVariable("49ShieldWidthSlab2",0.0015);
  // S
  // Control.addVariable("49ShieldWidthSlab2",0.56);
   // bare Mn
   // Control.addVariable("49ShieldWidthSlab2",0.015);
   //Mn+cd
    // Control.addVariable("49ShieldWidthSlab2",0.142);

   // 
  Control.addVariable("49ShieldMaterialSlab2","Void");  
  Control.addVariable("49ShieldWidthSlab3",5.08); 
  Control.addVariable("49ShieldMaterialSlab3","sbadMildSteel");  
  Control.addVariable("49ShieldWidthSlab4",0.70); 
  Control.addVariable("49ShieldMaterialSlab4","Void");  
  Control.addVariable("49ShieldWidthSlab5",5.08); 
  Control.addVariable("49ShieldMaterialSlab5","sbadMildSteel");  
  Control.addVariable("49ShieldWidthSlab6",0.70); 
  Control.addVariable("49ShieldMaterialSlab6","Void");  
  Control.addVariable("49ShieldWidthSlab7",5.08); 
  Control.addVariable("49ShieldMaterialSlab7","sbadMildSteel");  
  Control.addVariable("49ShieldWidthSlab8",0.70); 
  Control.addVariable("49ShieldMaterialSlab8","Void");  
  Control.addVariable("49ShieldWidthSlab9",3.00); 
  Control.addVariable("49ShieldMaterialSlab9","Stainless304");  
  Control.addVariable("49ShieldWidthSlab10",18.31); 
  Control.addVariable("49ShieldMaterialSlab10","H2O");  
  Control.addVariable("49ShieldWidthSlab11",3.00); 
  Control.addVariable("49ShieldMaterialSlab11","Stainless304");  
  Control.addVariable("49ShieldWidthSlab12",0.70); 
  Control.addVariable("49ShieldMaterialSlab12","Void");  
  Control.addVariable("49ShieldWidthSlab13",5.08); 
  Control.addVariable("49ShieldMaterialSlab13","sbadMildSteel");  
  Control.addVariable("49ShieldWidthSlab14",0.70); 
  Control.addVariable("49ShieldMaterialSlab14","Void");  
  Control.addVariable("49ShieldWidthSlab15",2.50); 
  Control.addVariable("49ShieldMaterialSlab15","sbadSulphurSteel");
  Control.addVariable("49ShieldWidthSlab16",19.80); 
  Control.addVariable("49ShieldMaterialSlab16","H2O");  
  Control.addVariable("49ShieldWidthSlab17",2.50); 
  Control.addVariable("49ShieldMaterialSlab17","Stainless304");  
  Control.addVariable("49ShieldWidthSlab18",0.70); 
  Control.addVariable("49ShieldMaterialSlab18","Void");  
  Control.addVariable("49ShieldWidthSlab19",5.08); 
  Control.addVariable("49ShieldMaterialSlab19","sbadMildSteel");  
  Control.addVariable("49ShieldWidthSlab20",0.70); 
  Control.addVariable("49ShieldMaterialSlab20","Void");  
  Control.addVariable("49ShieldWidthSlab21",5.08); 
  Control.addVariable("49ShieldMaterialSlab21","sbadMildSteel");  
  Control.addVariable("49ShieldWidthSlab22",30.70); 
  Control.addVariable("49ShieldMaterialSlab22","Void");  
  Control.addVariable("49ShieldWidthSlab23",15.24); 
  Control.addVariable("49ShieldMaterialSlab23","Concrete");  
  Control.addVariable("49ShieldWidthSlab24",2.00); 
  Control.addVariable("49ShieldMaterialSlab24","Void");  
  Control.addVariable("49ShieldWidthSlab25",15.24); 
  Control.addVariable("49ShieldMaterialSlab25","Concrete");  

  //Nestor side

  Control.addVariable("NestorXStep",91.45);  
  Control.addVariable("NestorYStep",-1.45);  
  Control.addVariable("NestorZStep",95.5);  

  Control.addVariable("NestorWidth",182.9);   
  Control.addVariable("NestorHeight",191.0);  

  Control.addVariable("NestorNSlab",6);
  Control.addVariable("NestorThick0",0.70); 
  Control.addVariable("NestorMat0","Void"); 
  Control.addVariable("NestorThick1",15.0); 
  Control.addVariable("NestorMat1","Graphite"); 
  Control.addVariable("NestorThick2",0.60); 
  Control.addVariable("NestorMat2","Void");  
  Control.addVariable("NestorThick3",5.08); 
  Control.addVariable("NestorMat3","sbadLead");  
  Control.addVariable("NestorThick4",0.70); 
  Control.addVariable("NestorMat4","Void");  
  Control.addVariable("NestorThick5",3.18); 
  Control.addVariable("NestorMat5","sbadMildSteel");  
  Control.addVariable("NestorAlWindowRadius",56.06); 


  // Fission Plate

  // Control.addVariable("49FissionPlateOffSetX",94.0);  
 Control.addVariable("49FissionPlateOffSetX",2.5);  

  Control.addVariable("49FissionPlateOffSetY",-1.45);   
 //Control.addVariable("49FissionPlateOffSetY",0.0);   

 // Control.addVariable("49FissionPlateOffSetZ",88.9);  
 Control.addVariable("49FissionPlateOffSetZ",-6.6);  


  // Control.addVariable("49FissionPlateLengthSlab",180.0);   
  // Control.addVariable("49FissionPlateHeightSlab",190.0);  
  // Control.addVariable("49FissionPlateLength",180.0);   
  // Control.addVariable("49FissionPlateHeight",190.0);  

 Control.addVariable("49FissionPlateLengthL",93.95);   
  Control.addVariable("49FissionPlateLengthR",88.95);   

  // correct but to avoid geom error it is the same as other slabs
  // Control.addVariable("35FissionPlateLengthSlab",182.9);   
   Control.addVariable("49FissionPlateHeightU",102.1);  
   Control.addVariable("49FissionPlateHeightD",88.9);  

  Control.addVariable("49FissionPlateTemperatureSlab",300.0); 

  // Control.addVariable("49FissionPlateFrontSlabN",7);
  // Control.addVariable("49FissionPlateWidthSlab0",1.2); 
  // Control.addVariable("49FissionPlateMaterialSlab0","sbadMildSteel");   
  // Control.addVariable("49FissionPlateWidthSlab1",0.1); 
  // Control.addVariable("49FissionPlateMaterialSlab1","Void");  
  // Control.addVariable("49FissionPlateWidthSlab2",0.1); 
  // Control.addVariable("49FissionPlateMaterialSlab2","sbadMildSteel");  
  // Control.addVariable("49FissionPlateWidthSlab3",0.1); 
  // Control.addVariable("49FissionPlateMaterialSlab3",4);  
  // Control.addVariable("49FissionPlateWidthSlab4",0.1); 
  // Control.addVariable("49FissionPlateMaterialSlab4",4);  
  // Control.addVariable("49FissionPlateWidthSlab5",0.1); 
  // Control.addVariable("49FissionPlateMaterialSlab5","sbadMildSteel");  
  // Control.addVariable("49FissionPlateWidthSlab6",1.2); 
  // Control.addVariable("49FissionPlateMaterialSlab6","sbadMildSteel");  



  Control.addVariable("49FissionPlateFrontSlabN",6);

  Control.addVariable("49FissionPlateLengthSlab0",119.0);   
  Control.addVariable("49FissionPlateHeightSlab0",142.0);  
  Control.addVariable("49FissionPlateWidthSlab0",1.2); 
  Control.addVariable("49FissionPlateMaterialSlab0","sbadMildSteel");   

  Control.addVariable("49FissionPlateLengthSlab1",119.0);   
  Control.addVariable("49FissionPlateHeightSlab1",142.0);  
  Control.addVariable("49FissionPlateWidthSlab1",0.1); 
  Control.addVariable("49FissionPlateMaterialSlab1","Void");  

  Control.addVariable("49FissionPlateLengthSlab2",119.0);   
  Control.addVariable("49FissionPlateHeightSlab2",142.0);  
  Control.addVariable("49FissionPlateWidthSlab2",0.1); 
  Control.addVariable("49FissionPlateMaterialSlab2","sbadMildSteel");

  Control.addVariable("49FissionPlateLengthSlab3",104.5);   
  Control.addVariable("49FissionPlateHeightSlab3",102.88);    
  Control.addVariable("49FissionPlateWidthSlab3",0.2); 
  Control.addVariable("49FissionPlateMaterialSlab3","Uranium");
  
  // Control.addVariable("49FissionPlateWidthSlab4",0.1); 
  // Control.addVariable("49FissionPlateMaterialSlab4",4);
 
  Control.addVariable("49FissionPlateLengthSlab4",119.0);   
  Control.addVariable("49FissionPlateHeightSlab4",142.0);   
  Control.addVariable("49FissionPlateWidthSlab4",0.1); 
  Control.addVariable("49FissionPlateMaterialSlab4","sbadMildSteel");  

  Control.addVariable("49FissionPlateLengthSlab5",119.0);   
  Control.addVariable("49FissionPlateHeightSlab5",142.0);  
  Control.addVariable("49FissionPlateWidthSlab5",1.2); 
  Control.addVariable("49FissionPlateMaterialSlab5","sbadMildSteel");  


  Control.addVariable("49Detector0Active1",1);   
  Control.addVariable("49DetectorOffSetX1",2.5);  
  Control.addVariable("49DetectorOffSetY1",5.1);  
  Control.addVariable("49DetectorOffSetZ1",-6.6);  

  Control.addVariable("49DetectorActive2",1);   
  Control.addVariable("49DetectorOffSetX2",2.5);  
  Control.addVariable("49DetectorOffSetY2",10.18001);  
  Control.addVariable("49DetectorOffSetZ2",-6.6);  

  Control.addVariable("49DetectorActive3",1);   
  Control.addVariable("49DetectorOffSetX3",2.5);  
  Control.addVariable("49DetectorOffSetY3",15.96);  
  Control.addVariable("49DetectorOffSetZ3",-6.6);  

  Control.addVariable("49DetectorActive4",1);   
  Control.addVariable("49DetectorOffSetX4",2.5);  
  Control.addVariable("49DetectorOffSetY4",21.74);  
  Control.addVariable("49DetectorOffSetZ4",-6.6);  

  Control.addVariable("49DetectorActive5",1);   
  Control.addVariable("49DetectorOffSetX5",2.5);  
  Control.addVariable("49DetectorOffSetY5",25.44);  
  Control.addVariable("49DetectorOffSetZ5",-6.6);  

  Control.addVariable("49DetectorActive6",0);   
  Control.addVariable("49DetectorOffSetX6",2.5);  
  Control.addVariable("49DetectorOffSetY6",27.44);  
  Control.addVariable("49DetectorOffSetZ6",-6.6);  

  Control.addVariable("49DetectorActive7",0);   
  Control.addVariable("49DetectorOffSetX7",2.5);  
  Control.addVariable("49DetectorOffSetY7",29.44);  
  Control.addVariable("49DetectorOffSetZ7",-6.6);  

  Control.addVariable("49DetectorActive8",0);   
  Control.addVariable("49DetectorOffSetX8",2.5);  
  Control.addVariable("49DetectorOffSetY8",29.94);  
  Control.addVariable("49DetectorOffSetZ8",-6.6);  

  Control.addVariable("49DetectorActive9",1);   
  Control.addVariable("49DetectorOffSetX9",2.5);  
  Control.addVariable("49DetectorOffSetY9",32.94);  
  Control.addVariable("49DetectorOffSetZ9",-6.6);  

  Control.addVariable("49DetectorActive10",0);   
  Control.addVariable("49DetectorOffSetX10",2.5);  
  Control.addVariable("49DetectorOffSetY10",34.44);  
  Control.addVariable("49DetectorOffSetZ10",-6.6);  

  Control.addVariable("49DetectorActive11",0);   
  Control.addVariable("49DetectorOffSetX11",2.5);  
  Control.addVariable("49DetectorOffSetY11",38.94);  
  Control.addVariable("49DetectorOffSetZ11",-6.6);  

  Control.addVariable("49DetectorActive12",1);   
  Control.addVariable("49DetectorOffSetX12",2.5);  
  Control.addVariable("49DetectorOffSetY12",40.44);  
  Control.addVariable("49DetectorOffSetZ12",-6.6);  

  Control.addVariable("49DetectorActive13",0);   
  Control.addVariable("49DetectorOffSetX13",2.5);  
  Control.addVariable("49DetectorOffSetY13",43.44);  
  Control.addVariable("49DetectorOffSetZ13",-6.6);  

  Control.addVariable("49DetectorActive14",1);   
  Control.addVariable("49DetectorOffSetX14",2.5);  
  Control.addVariable("49DetectorOffSetY14",46.75);  
  Control.addVariable("49DetectorOffSetZ14",-6.6);  

  Control.addVariable("49DetectorActive15",1);   
  Control.addVariable("49DetectorOffSetX15",2.5);  
  Control.addVariable("49DetectorOffSetY15",52.53);  
  Control.addVariable("49DetectorOffSetZ15",-6.6);  

  Control.addVariable("49DetectorActive16",1);   
  Control.addVariable("49DetectorOffSetX16",2.5);  
  Control.addVariable("49DetectorOffSetY16",55.73);  
  Control.addVariable("49DetectorOffSetZ16",-6.6);  

  Control.addVariable("49DetectorActive17",1);   
  Control.addVariable("49DetectorOffSetX17",2.5);  
  Control.addVariable("49DetectorOffSetY17",57.33);  
  Control.addVariable("49DetectorOffSetZ17",-6.6);  

  Control.addVariable("49DetectorActive18",1);   
  Control.addVariable("49DetectorOffSetX18",2.5);  
  Control.addVariable("49DetectorOffSetY18",59.73);  
  Control.addVariable("49DetectorOffSetZ18",-6.6);  

  Control.addVariable("49DetectorActive19",1);   
  Control.addVariable("49DetectorOffSetX19",2.5);  
  Control.addVariable("49DetectorOffSetY19",63.73);  
  Control.addVariable("49DetectorOffSetZ19",-6.6);  

  Control.addVariable("49DetectorActive20",1);   
  Control.addVariable("49DetectorOffSetX20",2.5);  
  Control.addVariable("49DetectorOffSetY20",65.83);  
  Control.addVariable("49DetectorOffSetZ20",-6.6);  

  Control.addVariable("49DetectorActive21",1);   
  Control.addVariable("49DetectorOffSetX21",2.5);  
  Control.addVariable("49DetectorOffSetY21",71.43);  
  Control.addVariable("49DetectorOffSetZ21",-6.6);  

  Control.addVariable("49DetectorActive22",1);   
  Control.addVariable("49DetectorOffSetX22",2.5);  
  Control.addVariable("49DetectorOffSetY22",74.53);  
  Control.addVariable("49DetectorOffSetZ22",-6.6);  



  Control.addVariable("49DetectorLayerN",1);

 // Rh
  Control.addVariable("49DetectorDiameterLayer0",1.27);   
  Control.addVariable("49DetectorWidthLayer0",0.0015); 
  Control.addVariable("49DetectorMaterialLayer0","Rh");   



  // S cast
  // Control.addVariable("49DetectorDiameterLayer0",5.1);   
  //  Control.addVariable("49DetectorWidthLayer0",0.56); 
  // Control.addVariable("49DetectorMaterialLayer0",101);   


  // Cd
  // Control.addVariable("49DetectorDiameterLayer0",1.27);   
  // Control.addVariable("49DetectorWidthLayer0",0.127); //0.0 
  // Control.addVariable("49DetectorMaterialLayer0",106);   

  // Mn
  // Control.addVariable("49DetectorDiameterLayer1",1.27);   
  // Control.addVariable("49DetectorWidthLayer1",0.015); 
  // Control.addVariable("49DetectorMaterialLayer1",104);   




 
  // S cast



 //Cd
  Control.addVariable("49DetectorDiameterLayer2",1.27);   
  Control.addVariable("49DetectorWidthLayer2",0.127); //0.0 
  Control.addVariable("49DetectorMaterialLayer2","Cadnium");   
  // Au
  Control.addVariable("49DetectorDiameterLayer4",1.27);   
  Control.addVariable("49DetectorWidthLayer4",0.0005); 
  Control.addVariable("49DetectorMaterialLayer4","Gold");   

  // S pressed
  Control.addVariable("49DetectorDiameterLayerX",3.81);   
  Control.addVariable("49DetectorWidthLayerX",0.241); 
  Control.addVariable("49DetectorMateriaLayerX","Rh");   


  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////      75        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



  Control.addVariable("75ShieldOffSetX",90.0);  
  Control.addVariable("75ShieldOffSetY",0.0);  
  Control.addVariable("75ShieldOffSetZ",90.0);  

  Control.addVariable("75ShieldTemperatureSlab",300.0); 

  Control.addVariable("75ShieldFrontSlabN",9);

  Control.addVariable("75ShieldWidthSlab0",1.91); 
  Control.addVariable("75ShieldLengthSlab0",180.);   
  Control.addVariable("75ShieldHeightSlab0",180.);  
  Control.addVariable("75ShieldMaterialSlab0","sbadLead");
 
  Control.addVariable("75ShieldWidthSlab1",12.1); 
  // fictitiously divide water layer
  Control.addVariable("75ShieldLengthSlab1",68.5);   
  Control.addVariable("75ShieldHeightSlab1",68.5);  
  // Control.addVariable("75ShieldLengthSlab1",180.);   
  // Control.addVariable("75ShieldHeightSlab1",180.);  
  Control.addVariable("75ShieldMaterialSlab1","Uranium"); 

  Control.addVariable("75ShieldWidthSlab2",5.9);
  Control.addVariable("75ShieldLengthSlab2",68.5);   
  Control.addVariable("75ShieldHeightSlab2",68.5);  
  Control.addVariable("75ShieldMaterialSlab2","sbadMildSteel");
  
  Control.addVariable("75ShieldWidthSlab3",12.7); 
  // fictitiously divide water layer
  Control.addVariable("75ShieldLengthSlab3",68.5);   
  Control.addVariable("75ShieldHeightSlab3",68.5);  
  // Control.addVariable("75ShieldLengthSlab3",180.);   
  // Control.addVariable("75ShieldHeightSlab3",180.);  
  Control.addVariable("75ShieldMaterialSlab3","Uranium");  

  Control.addVariable("75ShieldWidthSlab4",22.5); 
  Control.addVariable("75ShieldLengthSlab4",68.5);   
  Control.addVariable("75ShieldHeightSlab4",68.5);  
  Control.addVariable("75ShieldMaterialSlab4",6);  

  Control.addVariable("75ShieldWidthSlab5",0.3); 
  Control.addVariable("75ShieldLengthSlab5",180.);   
  Control.addVariable("75ShieldHeightSlab5",180.);  
  Control.addVariable("75ShieldMaterialSlab5","sbadLead");
  
  Control.addVariable("75ShieldWidthSlab6",29.58); 
  Control.addVariable("75ShieldLengthSlab6",59.4);   
  Control.addVariable("75ShieldHeightSlab6",59.4);  
  Control.addVariable("75ShieldMaterialSlab6","Void");  

  Control.addVariable("75ShieldWidthSlab7",0.3); 
  Control.addVariable("75ShieldLengthSlab7",60.);   
  Control.addVariable("75ShieldHeightSlab7",60.);  
  Control.addVariable("75ShieldMaterialSlab7","sbadLead");

  // inventato!!!!
  Control.addVariable("75ShieldWidthSlab8",50.0); 
  Control.addVariable("75ShieldLengthSlab8",180.);   
  Control.addVariable("75ShieldHeightSlab8",180.);  
  Control.addVariable("75ShieldMaterialSlab8","Uranium");  


  //Nestor side

  Control.addVariable("75NestorSideOffSetX",90.0);  
  Control.addVariable("75NestorSideOffSetY",0.0);  
  Control.addVariable("75NestorSideOffSetZ",90.0);  

  Control.addVariable("75NestorSideLengthSlab",180.);   
  Control.addVariable("75NestorSideHeightSlab",180.);  
  Control.addVariable("75NestorSideTemperatureSlab",300.0); 

  Control.addVariable("75NestorSideFrontSlabN",2);
  // attenzione: definire spessori a partire da origine verso l'indietro
  Control.addVariable("75NestorSideWidthSlab0",2.0); 
  Control.addVariable("75NestorSideMaterialSlab0","Void"); 
  Control.addVariable("75NestorSideWidthSlab1",15.0); 
  Control.addVariable("75NestorSideMaterialSlab1",1); 


  // Sinbad Cave

  //Nestor side

  Control.addVariable("75CaveOffSetX",0.0);  
  Control.addVariable("75CaveOffSetY",0.0);  
  Control.addVariable("75CaveOffSetZ",0.0);  

  //L: left, R: right, U: up, D: down, I: in (towards nestor), O: out
  Control.addVariable("75CaveLengthL",90);   
  Control.addVariable("75CaveLengthR",90.);   
  Control.addVariable("75CaveHeightU",90.);  
  Control.addVariable("75CaveHeightD",90.);  
  Control.addVariable("75CaveWidthI",17.0);  
  Control.addVariable("75CaveWidthO",135.29);
  // filled up with water  
  Control.addVariable("75CaveMaterial","Uranium"); 
  Control.addVariable("75CaveTemperature",300.0); 


  Control.addVariable("75CaveSlabN",4);

  // incremental thickness
  Control.addVariable("75CaveLengthLSlab0",2.5);   
  Control.addVariable("75CaveLengthRSlab0",2.5);
  // h inventate   
  Control.addVariable("75CaveHeightUSlab0",2.5);  
  Control.addVariable("75CaveHeightDSlab0",2.5);  
  Control.addVariable("75CaveWidthISlab0",3.21);  
  Control.addVariable("75CaveWidthOSlab0",0.0);  
  Control.addVariable("75CaveMaterialSlab0","Void"); 
  // the al window in slab 0
  Control.addVariable("75CaveAlWindowRadius",56.06); 
  //

  Control.addVariable("75CaveLengthLSlab1",0.0);   
  Control.addVariable("75CaveLengthRSlab1",0.0);   
  Control.addVariable("75CaveHeightUSlab1",0.0);  
  Control.addVariable("75CaveHeightDSlab1",0.0);  
  Control.addVariable("75CaveWidthISlab1",5.38);  
  Control.addVariable("75CaveWidthOSlab1",0.0);  
  Control.addVariable("75CaveMaterialSlab1","Void"); 

  Control.addVariable("75CaveLengthLSlab2",-30.5);   
  Control.addVariable("75CaveLengthRSlab2",-30.5);   
  Control.addVariable("75CaveHeightUSlab2",0.0);  
  Control.addVariable("75CaveHeightDSlab2",0.0);  
  Control.addVariable("75CaveWidthISlab2",0.65);  
  Control.addVariable("75CaveWidthOSlab2",0.0);  
  Control.addVariable("75CaveMaterialSlab2","sbadLead");

  Control.addVariable("75CaveLengthLSlab3",0.0);   
  Control.addVariable("75CaveLengthRSlab3",0.0);   
  Control.addVariable("75CaveHeightUSlab3",0.0);  
  Control.addVariable("75CaveHeightDSlab3",0.0);  
  Control.addVariable("75CaveWidthISlab3",28.91);  
  Control.addVariable("75CaveWidthOSlab3",0.0);  
  Control.addVariable("75CaveMaterialSlab3",1); 



  Control.addVariable("75CaveWidthSlab0",2.0); 

  Control.addVariable("75CaveMaterialSlab0","Void"); 


  Control.addVariable("75CaveWidthSlab1",15.0); 
  Control.addVariable("75CaveMaterialSlab1",1); 





  // Fission Plate

  Control.addVariable("75FissionPlateOffSetX",90.0);  
  Control.addVariable("75FissionPlateOffSetY",-1.40);  
  Control.addVariable("75FissionPlateOffSetZ",90.0);  

  // Control.addVariable("75FissionPlateLengthSlab",180.0);   
  // Control.addVariable("75FissionPlateHeightSlab",190.0);  
  Control.addVariable("75FissionPlateTemperatureSlab",300.0); 

  Control.addVariable("75FissionPlateFrontSlabN",8);

  Control.addVariable("75FissionPlateLengthSlab0",47.5);   
  Control.addVariable("75FissionPlateHeightSlab0",68.5);  
  Control.addVariable("75FissionPlateWidthSlab0",0.40); 
  Control.addVariable("75FissionPlateMaterialSlab0","sbadLead");   

  Control.addVariable("75FissionPlateLengthSlab1",40.2);   
  Control.addVariable("75FissionPlateHeightSlab1",63.5);  
  Control.addVariable("75FissionPlateWidthSlab1",0.1); 
  Control.addVariable("75FissionPlateMaterialSlab1","Void");  

  Control.addVariable("75FissionPlateLengthSlab2",40.2);   
  Control.addVariable("75FissionPlateHeightSlab2",63.5);  
  Control.addVariable("75FissionPlateWidthSlab2",0.4); 
  Control.addVariable("75FissionPlateMaterialSlab2",3);  

  Control.addVariable("75FissionPlateLengthSlab3",40.2);   
  Control.addVariable("75FissionPlateHeightSlab3",63.5);  
  Control.addVariable("75FissionPlateWidthSlab3",0.0); 
  Control.addVariable("75FissionPlateMaterialSlab3",3);  

  Control.addVariable("75FissionPlateLengthSlab4",40.2);   
  Control.addVariable("75FissionPlateHeightSlab4",63.5);  
  Control.addVariable("75FissionPlateWidthSlab4",0.0); 
  Control.addVariable("75FissionPlateMaterialSlab4",3);  

  Control.addVariable("75FissionPlateLengthSlab5",40.2);   
  Control.addVariable("75FissionPlateHeightSlab5",63.5);  
  Control.addVariable("75FissionPlateWidthSlab5",0.0); 
  Control.addVariable("75FissionPlateMaterialSlab5",3);  

  Control.addVariable("75FissionPlateLengthSlab6",40.2);   
  Control.addVariable("75FissionPlateHeightSlab6",63.5);  
  Control.addVariable("75FissionPlateWidthSlab6",0.1); 
  Control.addVariable("75FissionPlateMaterialSlab6","Void");  

  Control.addVariable("75FissionPlateLengthSlab7",47.5);   
  Control.addVariable("75FissionPlateHeightSlab7",68.5); 
  Control.addVariable("75FissionPlateWidthSlab7",0.4); 
  Control.addVariable("75FissionPlateMaterialSlab7","sbadLead");
  

  // Detectors

  Control.addVariable("75DetectorPositionN",10);  

  Control.addVariable("75DetectorActive0",1);   
  Control.addVariable("75DetectorOffSetX0",0.0);  
  Control.addVariable("75DetectorOffSetY0",1.91);  
  Control.addVariable("75DetectorOffSetZ0",0.0);  

  Control.addVariable("75DetectorActive1",1);   
  Control.addVariable("75DetectorOffSetX1",0.0);  
  Control.addVariable("75DetectorOffSetY1",7.41);  
  Control.addVariable("75DetectorOffSetZ1",0.0);  

  Control.addVariable("75DetectorActive2",1);   
  Control.addVariable("75DetectorOffSetX2",0.0);  
  Control.addVariable("75DetectorOffSetY2",12.41);  
  Control.addVariable("75DetectorOffSetZ2",0.0);  

  Control.addVariable("75DetectorActive3",1);   
  Control.addVariable("75DetectorOffSetX3",0.0);  
  Control.addVariable("75DetectorOffSetY3",14.01);  
  Control.addVariable("75DetectorOffSetZ3",0.0);  

  Control.addVariable("75DetectorActive4",1);   
  Control.addVariable("75DetectorOffSetX4",0.0);  
  Control.addVariable("75DetectorOffSetY4",19.91);  
  Control.addVariable("75DetectorOffSetZ4",0.0);  

  Control.addVariable("75DetectorActive5",1);   
  Control.addVariable("75DetectorOffSetX5",0.0);  
  Control.addVariable("75DetectorOffSetY5",25.41);  
  Control.addVariable("75DetectorOffSetZ5",0.0);  

  Control.addVariable("75DetectorActive6",1);   
  Control.addVariable("75DetectorOffSetX6",0.0);  
  Control.addVariable("75DetectorOffSetY6",30.41);  
  Control.addVariable("75DetectorOffSetZ6",0.0);  

  Control.addVariable("75DetectorActive7",1);   
  Control.addVariable("75DetectorOffSetX7",0.0);  
  Control.addVariable("75DetectorOffSetY7",39.01);  
  Control.addVariable("75DetectorOffSetZ7",0.0);  

  Control.addVariable("75DetectorActive8",1);   
  Control.addVariable("75DetectorOffSetX8",0.0);  
  Control.addVariable("75DetectorOffSetY8",49.61);  
  Control.addVariable("75DetectorOffSetZ8",0.0);  

  Control.addVariable("75DetectorActive9",1);   
  Control.addVariable("75DetectorOffSetX9",0.0);  
  Control.addVariable("75DetectorOffSetY9",58.61);  
  Control.addVariable("75DetectorOffSetZ9",0.0);  

  // spectrometers

 // Control.addVariable("75DetectorActive0",1);   
 //  Control.addVariable("75DetectorOffSetX0",0.0);  
 //  Control.addVariable("75DetectorOffSetY0",37.01);  
 //  Control.addVariable("75DetectorOffSetZ0",0.0);  

 //  Control.addVariable("75DetectorActive1",1);   
 //  Control.addVariable("75DetectorOffSetX1",0.0);  
 //  Control.addVariable("75DetectorOffSetY1",58.61);  
 //  Control.addVariable("75DetectorOffSetZ1",0.0);  

  // Mn 

  // Control.addVariable("75DetectorActive0","Void");   
  // Control.addVariable("75DetectorOffSetX0",0.0);  
  // Control.addVariable("75DetectorOffSetY0",1.91);  
  // Control.addVariable("75DetectorOffSetZ0",0.0);  

  // Control.addVariable("75DetectorActive1",0);   
  // Control.addVariable("75DetectorOffSetX1",0.0);  
  // Control.addVariable("75DetectorOffSetY1",2.41);  
  // Control.addVariable("75DetectorOffSetZ1",0.0);  

  // Control.addVariable("75DetectorActive2",0);   
  // Control.addVariable("75DetectorOffSetX2",0.0);  
  // Control.addVariable("75DetectorOffSetY2",3.41);  
  // Control.addVariable("75DetectorOffSetZ2",0.0);  

  // Control.addVariable("75DetectorActive3",0);   
  // Control.addVariable("75DetectorOffSetX3",0.0);  
  // Control.addVariable("75DetectorOffSetY3",4.41);  
  // Control.addVariable("75DetectorOffSetZ3",0.0);  

  // Control.addVariable("75DetectorActive4",0);   
  // Control.addVariable("75DetectorOffSetX4",0.0);  
  // Control.addVariable("75DetectorOffSetY4",5.41);  
  // Control.addVariable("75DetectorOffSetZ4",0.0);  

  // Control.addVariable("75DetectorActive5",0);   
  // Control.addVariable("75DetectorOffSetX5",0.0);  
  // Control.addVariable("75DetectorOffSetY5",6.41);  
  // Control.addVariable("75DetectorOffSetZ5",0.0);  

  // Control.addVariable("75DetectorActive6",0);   
  // Control.addVariable("75DetectorOffSetX6",0.0);  
  // Control.addVariable("75DetectorOffSetY6",7.51);  
  // Control.addVariable("75DetectorOffSetZ6",0.0);  

  // Control.addVariable("75DetectorActive7",0);   
  // Control.addVariable("75DetectorOffSetX7",0.0);  
  // Control.addVariable("75DetectorOffSetY7",8.41);  
  // Control.addVariable("75DetectorOffSetZ7",0.0);  

  // Control.addVariable("75DetectorActive8",0);   
  // Control.addVariable("75DetectorOffSetX8",0.0);  
  // Control.addVariable("75DetectorOffSetY8",9.41);  
  // Control.addVariable("75DetectorOffSetZ8",0.0);  

  // Control.addVariable("75DetectorActive9",0);   
  // Control.addVariable("75DetectorOffSetX9",0.0);  
  // Control.addVariable("75DetectorOffSetY9",10.51);  
  // Control.addVariable("75DetectorOffSetZ9",0.0);  


  // Control.addVariable("75DetectorActive10",0);   
  // Control.addVariable("75DetectorOffSetX10",0.0);  
  // Control.addVariable("75DetectorOffSetY10",11.41);  
  // Control.addVariable("75DetectorOffSetZ10",0.0);  

  // Control.addVariable("75DetectorActive11",0);   
  // Control.addVariable("75DetectorOffSetX11",0.0);  
  // Control.addVariable("75DetectorOffSetY11",13.51);  
  // Control.addVariable("75DetectorOffSetZ11",0.0);  

  // Control.addVariable("75DetectorActive12",0);   
  // Control.addVariable("75DetectorOffSetX12",0.0);  
  // Control.addVariable("75DetectorOffSetY12",14.01);  
  // Control.addVariable("75DetectorOffSetZ12",0.0);  

  // Control.addVariable("75DetectorActive13",0);   
  // Control.addVariable("75DetectorOffSetX13",0.0);  
  // Control.addVariable("75DetectorOffSetY13",19.91);  
  // Control.addVariable("75DetectorOffSetZ13",0.0);  

  // Control.addVariable("75DetectorActive14",0);   
  // Control.addVariable("75DetectorOffSetX14",0.0);  
  // Control.addVariable("75DetectorOffSetY14",20.41);  
  // Control.addVariable("75DetectorOffSetZ14",0.0);  

  // Control.addVariable("75DetectorActive15",0);   
  // Control.addVariable("75DetectorOffSetX15",0.0);  
  // Control.addVariable("75DetectorOffSetY15",21.41);  
  // Control.addVariable("75DetectorOffSetZ15",0.0);  

  // Control.addVariable("75DetectorActive16",0);   
  // Control.addVariable("75DetectorOffSetX16",0.0);  
  // Control.addVariable("75DetectorOffSetY16",22.11);  
  // Control.addVariable("75DetectorOffSetZ16",0.0);  

  // Control.addVariable("75DetectorActive17",0);   
  // Control.addVariable("75DetectorOffSetX17",0.0);  
  // Control.addVariable("75DetectorOffSetY17",22.41);  
  // Control.addVariable("75DetectorOffSetZ17",0.0);  

  // Control.addVariable("75DetectorActive18",0);   
  // Control.addVariable("75DetectorOffSetX18",0.0);  
  // Control.addVariable("75DetectorOffSetY18",23.41);  
  // Control.addVariable("75DetectorOffSetZ18",0.0);  

  // Control.addVariable("75DetectorActive19",0);   
  // Control.addVariable("75DetectorOffSetX19",0.0);  
  // Control.addVariable("75DetectorOffSetY19",24.11);  
  // Control.addVariable("75DetectorOffSetZ19",0.0);  


  // Control.addVariable("75DetectorActive20",0);   
  // Control.addVariable("75DetectorOffSetX20",0.0);  
  // Control.addVariable("75DetectorOffSetY20",24.41);  
  // Control.addVariable("75DetectorOffSetZ20",0.0);  

  // Control.addVariable("75DetectorActive21",0);   
  // Control.addVariable("75DetectorOffSetX21",0.0);  
  // Control.addVariable("75DetectorOffSetY21",26.11);  
  // Control.addVariable("75DetectorOffSetZ21",0.0);  

  // Control.addVariable("75DetectorActive22",0);   
  // Control.addVariable("75DetectorOffSetX22",0.0);  
  // Control.addVariable("75DetectorOffSetY22",26.41);  
  // Control.addVariable("75DetectorOffSetZ22",0.0);  

  // Control.addVariable("75DetectorActive23",0);   
  // Control.addVariable("75DetectorOffSetX23",0.0);  
  // Control.addVariable("75DetectorOffSetY23",27.41);  
  // Control.addVariable("75DetectorOffSetZ23",0.0);  

  // Control.addVariable("75DetectorActive24",0);   
  // Control.addVariable("75DetectorOffSetX24",0.0);  
  // Control.addVariable("75DetectorOffSetY24",28.11);  
  // Control.addVariable("75DetectorOffSetZ24",0.0);  

  // Control.addVariable("75DetectorActive25",0);   
  // Control.addVariable("75DetectorOffSetX25",0.0);  
  // Control.addVariable("75DetectorOffSetY25",29.11);  
  // Control.addVariable("75DetectorOffSetZ25",0.0);  

  // Control.addVariable("75DetectorActive26",0);   
  // Control.addVariable("75DetectorOffSetX26",0.0);  
  // Control.addVariable("75DetectorOffSetY26",30.11);  
  // Control.addVariable("75DetectorOffSetZ26",0.0);  

  // Control.addVariable("75DetectorActive27",0);   
  // Control.addVariable("75DetectorOffSetX27",0.0);  
  // Control.addVariable("75DetectorOffSetY27",31.11);  
  // Control.addVariable("75DetectorOffSetZ27",0.0);  

  // Control.addVariable("75DetectorActive28",0);   
  // Control.addVariable("75DetectorOffSetX28",0.0);  
  // Control.addVariable("75DetectorOffSetY28",32.11);  
  // Control.addVariable("75DetectorOffSetZ28",0.0);  

  // Control.addVariable("75DetectorActive29",0);   
  // Control.addVariable("75DetectorOffSetX29",0.0);  
  // Control.addVariable("75DetectorOffSetY29",32.61);  
  // Control.addVariable("75DetectorOffSetZ29",0.0);  



  Control.addVariable("75DetectorLayerN",1);

  //////////////////  choose det   ///////////////


 // Rh
   Control.addVariable("75DetectorDiameterLayer0",1.27);   
   Control.addVariable("75DetectorWidthLayer0",0.0015); 
   Control.addVariable("75DetectorMaterialLayer0",103);   




  // Spectrometer

  // Control.addVariable("75DetectorDiameterLayer0",-4);   
  // Control.addVariable("75DetectorWidthLayer0",-4.72); 
  // Control.addVariable("75DetectorMaterialLayer0","Void");   


  // Control.addVariable("75DetectorDiameterLayer1",-4);   
  // Control.addVariable("75DetectorWidthLayer1",-4.72); 
  // Control.addVariable("75DetectorMaterialLayer1","Void");   



  // Mn
  // Control.addVariable("75DetectorDiameterLayer0",1.27);   
  // Control.addVariable("75DetectorWidthLayer0",0.015); 
  // Control.addVariable("75DetectorMaterialLayer0",104);   


 // In
   // Control.addVariable("75DetectorDiameterLayer0",3.8);   
   // Control.addVariable("75DetectorWidthLayer0",0.163); 
   // Control.addVariable("75DetectorMaterialLayer0",107);   



  // S cast

  // Control.addVariable("75DetectorDiameterLayer0",5.1);   
  //  Control.addVariable("75DetectorWidthLayer0",0.56); 
  // Control.addVariable("75DetectorMaterialLayer0",101);   







 


 //Cd
  Control.addVariable("75DetectorDiameterLayer2",1.27);   
  Control.addVariable("75DetectorWidthLayer2",0.127); //0.0 
  Control.addVariable("75DetectorMaterialLayer2",106);   
  // Au
  Control.addVariable("75DetectorDiameterLayer4",1.27);   
  Control.addVariable("75DetectorWidthLayer4",0.0005); 
  Control.addVariable("75DetectorMaterialLayer4",105);   

  // S pressed
  Control.addVariable("75DetectorDiameterLayerX",3.81);   
  Control.addVariable("75DetectorWidthLayerX",0.241); 
  Control.addVariable("75DetectorMateriaLayerX",102);   











  return;
}

}  // NAMESPACE setVariable
