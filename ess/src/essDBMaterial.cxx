/********************************************************************* 
  CombLayer : MNCPX Input builder
 
 * File:   essBuild/essDBMaterial.cxx
*
 * Copyright (c) 2004-2013 by Stuart Ansell
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>. 
 *
 ****************************************************************************/
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cmath>
#include <set>
#include <map>
#include <list>
#include <vector>
#include <string>
#include <algorithm>
#include <functional>
#include <iterator>
#include <numeric>
#include <boost/bind.hpp>
#include <boost/tuple/tuple.hpp>

#include "Exception.h"
#include "FileReport.h"
#include "GTKreport.h"
#include "NameStack.h"
#include "RegMethod.h"
#include "OutputLog.h"
#include "support.h"
#include "BaseVisit.h"
#include "BaseModVisit.h"
#include "Element.h"
#include "Zaid.h"
#include "MXcards.h"
#include "Material.h"
#include "DBMaterial.h"
#include "essDBMaterial.h"

namespace ModelSupport
{

void addESSMaterial()
  /*!
     Initialize the database of materials
   */
{
  //  return;
  ELog::RegMethod RegA("essDBMaterial[F]","addESSMaterial");

  const std::string MLib="hlib=.70h pnlib=70u";
  ModelSupport::DBMaterial& MDB = ModelSupport::DBMaterial::Instance();

  MonteCarlo::Material MObj;


  // dummy materials for VTK
  MObj.setMaterial(1, "VTK1", " 1001.70c 1.0 ","", MLib);  MObj.setDensity(-1);  MDB.resetMaterial(MObj);
  MObj.setMaterial(2, "VTK2", " 1001.70c 1.0 ","", MLib);  MObj.setDensity(-1);  MDB.resetMaterial(MObj);
  MObj.setMaterial(3, "VTK3", " 1001.70c 1.0 ","", MLib);  MObj.setDensity(-1);  MDB.resetMaterial(MObj);
  MObj.setMaterial(4, "VTK4", " 1001.70c 1.0 ","", MLib);  MObj.setDensity(-1);  MDB.resetMaterial(MObj);
  MObj.setMaterial(5, "VTK5", " 1001.70c 1.0 ","", MLib);  MObj.setDensity(-1);  MDB.resetMaterial(MObj);
  MObj.setMaterial(6, "VTK6", " 1001.70c 1.0 ","", MLib);  MObj.setDensity(-1);  MDB.resetMaterial(MObj);
  MObj.setMaterial(7, "VTK7", " 1001.70c 1.0 ","", MLib);  MObj.setDensity(-1);  MDB.resetMaterial(MObj);
  MObj.setMaterial(8, "VTK8", " 1001.70c 1.0 ","", MLib);  MObj.setDensity(-1);  MDB.resetMaterial(MObj);
  MObj.setMaterial(9, "VTK9", " 1001.70c 1.0 ","", MLib);  MObj.setDensity(-1);  MDB.resetMaterial(MObj);
  MObj.setMaterial(10, "VTK10", " 1001.70c 1.0 ","", MLib);  MObj.setDensity(-1);  MDB.resetMaterial(MObj);
  MObj.setMaterial(11, "VTK11", " 1001.70c 1.0 ","", MLib);  MObj.setDensity(-1);  MDB.resetMaterial(MObj);
  // end dummy materials for VTK



  MObj.setMaterial(12, "M12",
		   " 1002.70c  0.333333333 "
		   " 1005.70c  0.666666666 ", // why this line changes density?
		   "dpara.60t dortho.60t", MLib);
  ELog::EM << "setMaterial 12: why density get changed with 1005.70c?" << ELog::endDiag;
  MObj.setDensity(-0.1624*2);  // kbat: I have to *2 due to the line above
  MObj.setMXitem(1002,70,'c',"h","01002 01002");
  MDB.resetMaterial(MObj);

  // ESS LIQUIDH2 at 20 K kbat
  MObj.setMaterial(1001, "M01001", " 1001.70c 1.000000000 ","hpara.10t", MLib);
  MObj.setDensity(-7.0e-2);
  MDB.resetMaterial(MObj);

  // ESS LIQUIDD2 at 20 K !!! check this again !!!
  MObj.setMaterial(1002, "M01002", " 1002.70c 1.000000000 ","dpara.10t", MLib);
  MObj.setDensity(-0.163);
  MDB.resetMaterial(MObj);

  // ESS Generic light water !!! check this again !!!
  MObj.setMaterial(1011, "M01011", " 01001.70c 0.66666667 "
		                   " 08016.70c 0.33333333 ", "lwtr.10t", MLib);
  MObj.setDensity(-1.0);
  MDB.resetMaterial(MObj);

  // ESS VIENNA STANDARD MEAN OCEAN WATER kbat
  MObj.setMaterial(1015, "M01015",
		   " 01001.70c  0.666562842 "
		   " 01002.70c  0.000103824 "
		   " 08016.70c  0.332540192 "
		   " 08017.70c  0.000126332 "
		   " 08018.70c   0.000666810", "lwtr.10t", MLib);
  MObj.setMXitem(8018,70,'c',"n","model");
  // If there is no library for the isotope, one should put it without extension and then say "model" in the mx card:
  // but CombLayer does not allow rectords without extensions, therefore I am using 08018.80c here.
  //  std::cerr << "WARNING M01015: Alan's definition specifies 08018 instead of 08018.70c" << std::endl;
  MObj.setDensity(-1.0);
  MDB.resetMaterial(MObj);

  // ESS GENERIC HEAVY WATER !!! check this again !!!
  MObj.setMaterial(1021, "M01021",
		   " 01002.70c  0.666666667 "
		   " 08016.70c  0.333333333 ", "hwtr.10t", MLib);
  //  MObj.setMXitem(8018,70,'c',"n","model");
  MObj.setDensity(-1.1);
  MDB.resetMaterial(MObj);

  MObj.setMaterial(1995, "M01995",
		   " 01001.70c 0.995 "
		   " 01004.70c 0.005 ", "hpara.10t hortho.10t", MLib);
  MObj.setDensity(-0.07);
  MDB.resetMaterial(MObj);


  // "Real life" para-ortho mixture
  // 99.8% para, 0.2% ortho - approx ground state at 20 K
  // Temperature: 20 K
  // Reference: Self-calculated by KB based on the Alan's email 28 Nov 2013.  Works with
  // S(a,b) libraries on the cluster.
  MObj.setMaterial(1998, "M01998",
		   " 01001.70c 0.998 "
		   " 01004.70c 0.002 ", "hpara.10t hortho.10t", MLib);
  MObj.setDensity(-0.07);
  MDB.resetMaterial(MObj);
  

  // ESS He in monolith !!! check this again !!!
  MObj.setMaterial(2000, "M02000",
		   " 02003.70c  0.000001340 "
		   " 02004.70c  0.999998660 ", "", MLib);
  //  MObj.setMXitem(8018,70,'c',"n","model");
  MObj.setDensity(-1.74E-4);
  MDB.resetMaterial(MObj);

  // ESS He in PBW !!! check this again !!!
  MObj.setMaterial(2001, "M02001",
		   " 02003.70c  0.000001340 "
		   " 02004.70c  0.999998660 ", "", MLib);
  //  MObj.setMXitem(8018,70,'c',"n","model");
  MObj.setDensity(-1.72E-3);
  MDB.resetMaterial(MObj);

  // ESS He in Target kbat
  MObj.setMaterial(2002, "M02002",
		   " 02003.71c  0.000001340 "
		   " 02004.71c  0.999998660 ", "", MLib);
  MObj.setDensity(-2.58E-4);
  MDB.resetMaterial(MObj);

  // ESS M04000 Be 300 K kbat
  MObj.setMaterial(4000,"M04000", "4009.70c  1.0 ", "be.10t",MLib);
  MObj.setDensity(-1.85);
  MObj.setMXitem(4009,70, 'c', "h", "model");
  MDB.resetMaterial(MObj); 

  // ESS M04001 Be !!! check this again !!!
  MObj.setMaterial(4001, "M04001", "4009.70c  1.0 ", "beme.31t",MLib);
  MObj.setDensity(-1.85);
  //  MObj.setMXitem(4009,70,'c',"h","model");
  MDB.resetMaterial(MObj); 

  // Be + 10% H2O coolant - from the Alan's email 13 Jun 2014 - checked by KB
  MObj.setMaterial(4010, "M04010",
		   " 1001.70c 0.055140611 "
		   " 4009.70c 0.91728908 "
		   " 8016.70c 0.027570305 ", "lwtr.10t be.10t", MLib);
  MObj.setMXitem(4009,70,'c',"h","model");
  MObj.setDensity(-1.76);
  MDB.resetMaterial(MObj);

  // Be + 10% D2O coolant - from the Alan's email 13 Jun 2014 - checked by KB
  MObj.setMaterial(4020, "M04020",
		   " 1002.70c 0.054605886 "
		   " 4009.70c 0.918091171 "
		   " 8016.70c 0.027302943 ", "hwtr.10t be.10t", MLib);
  MObj.setMXitem(4009,70,'c',"h","model");
  MObj.setDensity(-1.77);
  MDB.resetMaterial(MObj);


  // ESS M04011 Beryllane !!! check this again !!!
  MObj.setMaterial(4011, "M04011", 
		   " 01001.70c 0.666666667 "
		   " 04009.70c 0.333333333 ", "", MLib);
  MObj.setDensity(-0.65);
  //  MObj.setMXitem(4009,70,'c',"h","model");
  MDB.resetMaterial(MObj); 

  // Graphite
  MObj.setMaterial(6000, "M06000", "06000.70c 1.0", "grph.10t", MLib);
  MObj.setDensity(-2.3);
  MDB.resetMaterial(MObj);

  // ESS M08011 CONCRETEMINE !!! check this again !!!
  MObj.setMaterial(8011, "M08011", 
		   " 01001.70c 0.15 "
		   " 06000.70c 0.1  "
		   " 08016.70c 0.55 "
		   " 14028.70c 0.1  "
		   " 20040.70c 0.1 ", "lwtr.10t", MLib);
  MObj.setDensity(-2.3);
  MObj.setMXitem(6000,70,'c', "h", "06012");
  MDB.resetMaterial(MObj); 


  // ESS M13000 Al !!! check this again !!!
  MObj.setMaterial(13000, "M13000", " 13027.70c 1.0 ", "al27.12t", MLib);
  MObj.setDensity(-2.7);
  MDB.resetMaterial(MObj); 

  // ESS M13001 Al !!! check this again !!!
  MObj.setMaterial(13001, "M13001", " 13027.70c 1.0 ", "al27.10t", MLib);
  MObj.setDensity(-2.73);
  MDB.resetMaterial(MObj); 

  // ESS M13060 #81 Al-6061-T6 - same as 13061 but at normal T (and slightly different density) !!! isotope composition checked, check other numbers again !!!
  MObj.setMaterial(13060,"M13060"," 12024.70c  0.008812087 "
		   "12025.70c 0.001115595 12026.70c  0.001228270 "
		   "13027.70c 0.977848644 14028.70c  0.005342094 "
		   "14029.70c 0.000271383 14030.70c  0.000179107 "
		   "22046.70c 0.000035050 22047.70c  0.000031608 "
		   "22048.70c  0.000313196 22049.70c  0.000022984 "
		   "22050.70c  0.000022007 24050.70c  0.000044183 "
		   "24052.70c  0.000852028 24053.70c  0.000096613 "
		   "24054.70c  0.000024049 25055.70c  0.000370161 "
		   "26054.70c  0.000099328 26056.70c  0.001559232 "
		   "26057.70c  0.000036009 26058.70c  0.000004792 "
		   "29063.70c  0.000811409 29065.70c  0.000361995 "
		   "30000.70c  0.000518176 ","al27.12t",MLib);
  MObj.setDensity(-2.70);
  MDB.resetMaterial(MObj);

  // ESS M13061  Al-6061-T6 - same as 13060, but at 20 K (and slightly different density) !!! isotope composition checked, check other numbers again !!!
  MObj.setMaterial(13061,"M13061"," 12024.70c  0.008812087 "
		   "12025.70c 0.001115595 12026.70c  0.001228270 "
		   "13027.70c 0.977848644 14028.70c  0.005342094 "
		   "14029.70c 0.000271383 14030.70c  0.000179107 "
		   "22046.70c 0.000035050 22047.70c  0.000031608 "
		   "22048.70c  0.000313196 22049.70c  0.000022984 "
		   "22050.70c  0.000022007 24050.70c  0.000044183 "
		   "24052.70c  0.000852028 24053.70c  0.000096613 "
		   "24054.70c  0.000024049 25055.70c  0.000370161 "
		   "26054.70c  0.000099328 26056.70c  0.001559232 "
		   "26057.70c  0.000036009 26058.70c  0.000004792 "
		   "29063.70c  0.000811409 29065.70c  0.000361995 "
		   "30000.70c  0.000518176 ","al27.10t",MLib);
  MObj.setDensity(-2.73);
  MDB.resetMaterial(MObj);


  // ESS Iron - as in the Alan's table
  MObj.setMaterial(26000, "M26000",
		   " 26054.70c  0.058450000 "
		   " 26056.70c  0.917540000 "
		   " 26057.70c  0.021190000 "
		   " 26058.70c  0.002820000 ", "fe56.12t", MLib);
  //  MObj.setMXitem(6000, 70, 'c', "h", "06012");
  MObj.setDensity(-7.85);
  MDB.resetMaterial(MObj);

  // ESS STANDARD GRADE CARBON STEEL AISI 1005 kbat
  MObj.setMaterial(26005, "M26005",
		   " 06000.70c  0.002781585 "
		   " 15031.70c  0.000719079 "
		   " 16032.70c  0.000824765 "
		   " 16033.70c  0.000006512 "
		   " 16034.70c  0.000036901 "
		   " 16036.70c  0.000000087 "
		   " 25055.70c  0.003547362 "
		   " 26054.70c  0.057987293 "
		   " 26056.70c  0.910276486 "
		   " 26057.70c  0.021022254 "
		   " 26058.70c  0.002797676 ", "fe56.12t", MLib);
  MObj.setMXitem(6000, 70, 'c', "h", "06012");
  MObj.setDensity(-7.85);
  MDB.resetMaterial(MObj);

  // FE + 10%H2O - from the Alan's email 16 Jun 2014
  MObj.setMaterial(26010, "M26010",
		   " 01001.70c 0.077534884 "
		   " 08016.70c 0.038767442 "
		   " 26054.70c 0.051652129 "
		   " 26056.70c 0.810827964 "
		   " 26057.70c 0.018725554 "
		   " 26058.70c 0.002492027 ", "lwtr.10t fe56.12t", MLib);
  MObj.setDensity(-7.17);
  MDB.resetMaterial(MObj);

  // FE + 10%D2O - from the Alan's email 16 Jun 2014
  MObj.setMaterial(26020, "M26020",
		   " 01002.70c 0.076810269 "
		   " 08016.70c 0.038405135 "
		   " 26054.70c 0.051715660 "
		   " 26056.70c 0.811825257 "
		   " 26057.70c 0.018748586 "
		   " 26058.70c 0.002495093 ", "hwtr.10t fe56.12t", MLib);
  MObj.setDensity(-7.18);
  MDB.resetMaterial(MObj);

  // ESS SS316L kbat
  MObj.setMaterial(26316, "M26316",
		   " 06000.70c  0.001392603 "
		   " 14028.70c  0.007323064 "
		   " 14029.70c  0.000372017 "
		   " 14030.70c  0.000245523 "
		   " 15031.70c  0.000360008 "
		   " 16032.70c  0.000165168 "
		   " 16033.70c  0.000001304 "
		   " 16034.70c  0.000007390 "
		   " 16036.70c  0.000000017 "
		   " 24050.70c  0.007920331 "
		   " 24052.70c  0.152735704 "
		   " 24053.70c  0.017319003 "
		   " 24054.70c  0.004311066 "
		   " 25055.70c  0.018267327 "
		   " 26054.70c  0.038344779 "
		   " 26056.70c  0.601931034 "
		   " 26057.70c  0.013901213 "
		   " 26058.70c  0.001849996 "
		   " 27059.70c  0.000283816 "
		   " 28058.70c  0.080834464 "
		   " 28060.70c  0.031137291 "
		   " 28061.70c  0.001353516 "
		   " 28062.70c  0.004315603 "
		   " 28064.70c  0.001099057 "
		   " 42092.70c  0.002145890 "
		   " 42094.70c  0.001341000 "
		   " 42095.70c  0.002310064 "
		   " 42096.70c  0.002423388 "
		   " 42097.70c  0.001388944 "
		   " 42098.70c  0.003514494 "
		   " 42100.70c  0.001404926 ", "fe56.12t", MLib);
  MObj.setMXitem(6000, 70, 'c', "h", "06012");
  MObj.setDensity(-7.85);
  MDB.resetMaterial(MObj);


  // ESS  SS316L kbat
  MObj.setMaterial(26317, "M26317",
		   " 06000.71c  0.001392603 "
		   " 14028.71c  0.007323064 "
		   " 14029.71c  0.000372017 "
		   " 14030.71c  0.000245523 "
		   " 15031.71c  0.000360008 "
		   " 16032.71c  0.000165168 "
		   " 16033.71c  0.000001304 "
		   " 16034.71c  0.000007390 "
		   " 16036.71c  0.000000017 "
		   " 24050.71c  0.007920331 "
		   " 24052.71c  0.152735704 "
		   " 24053.71c  0.017319003 "
		   " 24054.71c  0.004311066 "
		   " 25055.71c  0.018267327 "
		   " 26054.71c  0.038344779 "
		   " 26056.71c  0.601931034 "
		   " 26057.71c  0.013901213 "
		   " 26058.71c  0.001849996 "
		   " 27059.71c  0.000283816 "
		   " 28058.71c  0.080834464 "
		   " 28060.71c  0.031137291 "
		   " 28061.71c  0.001353516 "
		   " 28062.71c  0.004315603 "
		   " 28064.71c  0.001099057 "
		   " 42092.71c  0.002145890 "
		   " 42094.71c  0.001341000 "
		   " 42095.71c  0.002310064 "
		   " 42096.71c  0.002423388 "
		   " 42097.71c  0.001388944 "
		   " 42098.71c  0.003514494 "
		   " 42100.71c  0.001404926 ", "fe56.14t", MLib);
  MObj.setMXitem(6000, 71, 'c', "h", "06012");
  MObj.setDensity(-7.76);
  MDB.resetMaterial(MObj);

  // Nickel from the Alan's table
  MObj.setMaterial(28000, "M28000", 
		   " 28058.70c 0.680769000 "
		   " 28060.70c 0.262231000 "
		   " 28061.70c 0.011399000 "
		   " 28062.70c 0.036345000 "
		   " 28064.70c 0.009256000 ",
		   "", MLib);
  MObj.setDensity(-8.90);
  MDB.resetMaterial(MObj);

  // Zn - kbat - not in the Alan's table 
  MObj.setMaterial(30000, "M30000",
		   " 30064.70c 0.486 "
		   " 30066.70c 0.279 "
		   " 30067.70c 0.041 "
		   " 30068.70c 0.188 "
		   " 30070.70c 0.06  ", "", MLib);
  MObj.setDensity(-7.14); // wiki
  MDB.resetMaterial(MObj);
  ELog::EM << "!!! essDBMaterial: check extensions for Zn !!!" << ELog::endDiag;

  // Zirconium - Alan - from webelements.
  // Source: M. Nomura, K. Kogure, M. Okamoto, Int. J. Mass Spectrom. Ion Phys. 50, 219 (1983)
  MObj.setMaterial(40000, "M40000",
		   " 40090.70c 0.5145 "
		   " 40091.70c 0.1122 "
		   " 40092.70c 0.1715 "
		   " 40094.70c 0.1738 "
		   " 40096.70c 0.0280 ", "", MLib);
  MObj.setDensity(-6.49);
  MDB.resetMaterial(MObj);
  ELog::EM << "!!! essDBMaterial: Alan's density of Zr is not the same as in WebElements !!!" << ELog::endDiag;

  // Zr+10% vol H2O
  // 
  // Temperature = 300 K
  // Reference: calculated with mc-tools based on materials 01011 and 40000
  MObj.setMaterial(40010, "M40010",
		   " 01001.70c 0.137600519836 "
		   " 08016.70c 0.068800259918 "
		   " 40090.70c 0.408306798817 "
		   " 40091.70c 0.0890418325116 "
		   " 40092.70c 0.136102266272 "
		   " 40094.70c 0.137927544479 "
		   " 40096.70c 0.0222207781669 ", "lwtr.10t", MLib);
  MObj.setDensity(-5.941);
  MDB.resetMaterial(MObj);

  //M 0115 #76
  MObj.setMaterial(76,"EssH2O","1001.70c  0.666562842 "
		   " 1002.70c  0.000103824 8016.70c  0.332540192 "
		   " 8017.70c  0.000126332 8018.70c  0.000666810 ",
		   "",MLib);
  MObj.setDensity(-1.0);
  MDB.resetMaterial(MObj);
 
  // M0200 #77
  MObj.setMaterial(77,"EssHe"," 2003.70c  0.00000134 "
		   " 2004.70c  0.99999866 " ,
		   "",MLib);
  MObj.setDensity(-1.74e-4);
  MDB.resetMaterial(MObj);
  
  
  // // M1361 #82
  // MObj.setMaterial(82,"M1361",
  // 		   "12024.70c 0.008812087 12025.70c 0.001115595 "
  // 		   "12026.70c 0.001228270 13027.70c 0.977848644 "
  // 		   "14028.70c 0.005342094 14029.70c 0.000271383 "
  // 		   "14030.70c 0.000179107 2046.70c 0.000035050 "
  // 		   "22047.70c 0.000031608 2048.70c 0.000313196 "
  // 		   "22049.70c 0.000022984 2050.70c 0.000022007 "
  // 		   "24050.70c 0.000044183 4052.70c 0.000852028 "
  // 		   "24053.70c 0.000096613 4054.70c 0.000024049 "
  // 		   "25055.70c 0.000370161 6054.70c 0.000099328 "
  // 		   "26056.70c 0.001559232 6057.70c 0.000036009 "
  // 		   "26058.70c 0.000004792 9063.70c 0.000811409 "
  // 		   "29065.70c 0.000361995 0000.70c 0.000518176 ",
  // 		   "",MLib);
  // MObj.setDensity(-2.73);
  // MDB.resetMaterial(MObj);
  // return;

  // M2600 #83 
  MObj.setMaterial(83,"M2600","26054.70c  0.058450000 "
		   "26056.70c 0.917540000 26057.70c 0.021190000 "
		   "26058.70c 0.002820000 ","",MLib);
  MObj.setDensity(-7.85);
  MDB.resetMaterial(MObj);

  //M2616 #84   
  MObj.setMaterial(84,"M2616",
		   "06000.71c  0.001392603 14028.71c  0.007323064 "
		   "14029.71c  0.000372017 14030.71c  0.000245523 "
		   "15031.71c  0.000360008 16032.71c  0.000165168 "
		   "16033.71c  0.000001304 16034.71c  0.000007390 "
		   "16036.71c  0.000000017 24050.71c  0.007920331 "
		   "24052.71c  0.152735704 24053.71c  0.017319003 "
		   "24054.71c  0.004311066 25055.71c  0.018267327 "
		   "26054.71c  0.038344779 26056.71c  0.601931034 " 
		   "26057.71c  0.013901213 26058.71c  0.001849996 "
		   "27059.71c  0.000283816 28058.71c  0.080834464 "
		   "28060.71c  0.031137291 28061.71c  0.001353516 "
		   "28062.71c  0.004315603 28064.71c  0.001099057 ",
		   "",MLib);
  MObj.setDensity(-7.76);
  MDB.resetMaterial(MObj);  


  //ESS M74000 kbat - same ase 74001 but at 300 K
  // T = 300 K
  MObj.setMaterial(74000,"M74000",
		   "74180.50c  0.001200000 "
		   "74182.71c  0.265000000 "
		   "74183.71c  0.143100000 "
		   "74184.71c  0.306400000 "
		   "74186.71c  0.284300000 ",
		   "",MLib);
  MObj.setDensity(-19.3);
  MDB.resetMaterial(MObj);  

  //ESS M74001 kbat
  // T = 600 K
  MObj.setMaterial(74001,"M74001",
		   "74180.50c  0.001200000 "
		   "74182.71c  0.265000000 "
		   "74183.71c  0.143100000 "
		   "74184.71c  0.306400000 "
		   "74186.71c  0.284300000 ",
		   "",MLib);
  MObj.setDensity(-19.20);
  MDB.resetMaterial(MObj);  

  //ESS M82000 kbat  Natural lead rho=11.4 g/cm3w
  MObj.setMaterial(82000,"M82000", 
		   " 82204.70c 0.014       "
		   " 82206.70c 0.241000000 "
		   " 82207.70c 0.221000000 "
		   " 82208.70c 0.524000000 ",
		   " pb206.20t pb207.20t pb208.20t ",MLib);
  MObj.setDensity(-11.4);
  MDB.resetMaterial(MObj);  

  // Lead + 10% vol H2O
  // 
  // Temperature = 300 K
  // Reference: calculated with mc-tools based on materials 01011 and 82000
  MObj.setMaterial(82010, "M82010",
		   " 01001.70c 0.167788489053 "
		   " 08016.70c 0.0838942445267 "
		   " 82204.70c 0.0104764417299 "
		   " 82206.70c 0.180344461207 "
		   " 82207.70c 0.165378115879 "
		   " 82208.70c 0.392118247604 ", "lwtr.10t pb206.20t pb207.20t pb208.20t", MLib);
  MObj.setDensity(-10.36);
  MDB.resetMaterial(MObj);

  // Lead + 10% vol D2O
  // 
  // Temperature: 300 K
  // Reference: Calculated with mc-tools
  MObj.setMaterial(82020, "M82020",
		   " 01002.70c 0.167256705756 "
		   " 08016.70c 0.083628352878 "
		   " 82204.70c 0.0104876091791 "
		   " 82206.70c 0.180536700869 "
		   " 82207.70c 0.165554402042 "
		   " 82208.70c 0.392536229276 ", "hwtr.10t pb206.20t pb207.20t pb208.20t", MLib);
  MObj.setDensity(-10.371);
  MDB.resetMaterial(MObj);

  // Bismuth - as in the Alan's table
  MObj.setMaterial(83000, "M83000",
		   "83209.70c 1.0", "bism.20t", MLib);
  MObj.setDensity(-9.8);
  MDB.resetMaterial(MObj);  
 

  
  return;
}

} // NAMESPACE ModelSupport
