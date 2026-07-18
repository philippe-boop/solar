/*-------------------------------------------------------------------------------*/
/*  SOLAR - The solar thermal power plant simulator                              */
/*  https://github.com/bbopt/solar                                               */
/*                                                                               */
/*  Miguel Diago, Sebastien Le Digabel, Mathieu Lemyre-Garneau, Bastien Talgorn  */
/*                                                                               */
/*  Polytechnique Montreal / GERAD                                               */
/*  sebastien.le-digabel@polymtl.ca                                              */
/*                                                                               */
/*  This program is free software: you can redistribute it and/or modify it      */
/*  under the terms of the GNU Lesser General Public License as published by     */
/*  the Free Software Foundation, either version 3 of the License, or (at your   */
/*  option) any later version.                                                   */
/*                                                                               */
/*  This program is distributed in the hope that it will be useful, but WITHOUT  */
/*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        */
/*  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  */
/*  for more details.                                                            */
/*                                                                               */
/*  You should have received a copy of the GNU Lesser General Public License     */
/*  along with this program. If not, see <http://www.gnu.org/licenses/>.         */
/*                                                                               */
/*-------------------------------------------------------------------------------*/
#include "MoltenSalt.hpp"

//Coefficients for viscosity determination
//Viscosity equation :
//eta = a + b*T + c*T^2 + d*T^3
//a, b, c and d vary as a function of the molar concentration of NaNO3 and KNO3
//here x is the concentration of NaNO3.
//a = 0.0431x2 - 7.9905x + 382.68
#define a_A 0.0431
#define a_B (-7.9905)
#define a_C 382.68
//b x 10 ^2 = -0.019x2 + 3.5305x - 164.19
#define b_A (-0.019)
#define b_B 3.5305
#define b_C (-164.19)
//c x 10^5 = 0.028x2 - 5.2045x + 236.84
#define c_A 0.028
#define c_B (-5.2045)
#define c_C 236.84
//d x 10^8 = -0.0138x2 + 2.5561x - 114.23
#define d_A (-0.0138)
#define d_B 2.5561
#define d_C (-114.23)

/*------------------------------------------------------------------------------------*/
/*                                  constructor #1                                    */
/*P.B. 2026-06 : Modified this constructor to allow different types of salt to be used*/
/*------------------------------------------------------------------------------------*/
MoltenSalt::MoltenSalt ( double temp, double pres, double masf, int id ) :
  _temperature(temp), _pressure(pres), _massFlow(masf), _id(id) {
    
    switch (id) {
      
      //60% NaNO3 - 40% KNO3 (Solar salt). Used in instances 1-11.
      case 0:
      _meltingPoint     = 495.0;
      _maxOperatingTemp = 0.0;  //Max. operating temperature is the upper bound of receiver outlet temp.
      _costPerKg        = 1.2; //Old data
      _heatCapacity     = 1530.0;
      break;

      //60% NaNO3 - 40% KNO3 (Solar salt)
      case 1:
      _meltingPoint     = 495.0;
      _maxOperatingTemp = 565.0 + 273.15;
      _costPerKg        = 0.80;
      _heatCapacity     = 1530.0;
      break;
      
    //32.12wt%Li2CO3 – 33.36wt%Na2CO3 –34.52wt% K2CO3 (denoted as (LiNaK)2CO3)
      case 2:
      _meltingPoint     = 396.8 + 273.15;
      _maxOperatingTemp = 800.0 + 273.15;
      _costPerKg        = 2.50;    
      _heatCapacity     = 1610.0;
      break;
    
      //CuCl-KCl-NaCl without additive
      case 3:
      _meltingPoint     = 144.0 + 273.15;
      _maxOperatingTemp = 653.0 + 273.15;
      _costPerKg        = 0.8; //Not found yet
      _heatCapacity     = 740.0;                
      break;
    
      //CuCl-KCl-NaCl with 7% CaCl2 additive
      case 4:
      _meltingPoint     = 138.3 + 273.15;
      _maxOperatingTemp = 700.0 + 273.15;
      _costPerKg        = 0.8; //Not found yet          
      _heatCapacity     = 850.0;
      break;

      //KNO3–NaNO3–LiNO3–Ca(NO3)2*4H2O 6:1:2:2
      case 5:
      _meltingPoint     = 88.6  + 273.15;
      _maxOperatingTemp = 612.0 + 273.15;
      _costPerKg        = 1.0; //Not found yet
      _heatCapacity     = 1571.0; //Evaluated at 350.3 Celsius
      break;
      
      default: throw std::invalid_argument("Invalid molten salt ID");
    }

}

/*-------------------------------------------------------------*/

double MoltenSalt::computeViscosity ( double T ) {
/*-------------------------------------------------------------*/
  double R = 8.314; // J/(mol K)
  
  switch (_id) {
    
    //60% NaNO3 - 40% KNO3 (Solar salt)    
    case 0:
    case 1: {
      double C = PERCENT_MASS_NANO3;

      // calculate mol percent from percent mass
      double N_NaNO3 = C / MOL_MASS_NANO3;
      double N_KNO3  = (1 - C) / MOL_MASS_KNO3;
      double N_total = N_NaNO3 + N_KNO3;
      double C_Mol   = (N_NaNO3 / N_total)*100.0;

      double a = a_A * pow(C_Mol, 2.) + a_B * C_Mol + a_C;
      double b = b_A * pow(C_Mol, 2.) + b_B * C_Mol + b_C;
      double c = c_A * pow(C_Mol, 2.) + c_B * C_Mol + c_C;
      double d = d_A * pow(C_Mol, 2.) + d_B * C_Mol + d_C;
      return (a + b*T*1e-2 + c*pow(T, 2.0)*1e-5 + d*pow(T, 3.0)*1e-8) / 1000.0;}
    
    //32.12wt%Li2CO3 – 33.36wt%Na2CO3 –34.52wt% K2CO3 (denoted as (LiNaK)2CO3)
    case 2: return 0.0852 * std::exp( (3.51e4) / (R*T) ) /1000.0;

    //CuCl-KCl-NaCl without additive    
    case 3: return (20.763 - 0.0255*(T - 273.15)) / 1000.0;
    
    //CuCl-KCl-NaCl with 7% CaCl2 additive    
    case 4: return (23.228 - 0.0306*(T - 273.15)) / 1000.0;

    //KNO3–NaNO3–LiNO3–Ca(NO3)2*4H2O 6:1:2:2
    case 5: return 0.6492 * std::exp( (7450.0) / (R*T)) / 1000.0;
    
    default: throw std::invalid_argument("Invalid molten salt ID");
  }
}

double MoltenSalt::computeDensity      ( double T) {
  switch (_id) {
    
    //Solar Salt
    case 0: return 1840.0;
    case 1: return 1840.0;

    // (LiNaK)2CO3
    case 2: return 1000.0*(2.27 - 4.34e-4 * (T - 273.15));

    //CuCl-KCl-NaCl without additive
    case 3: return 1000.0 * (3.1175 - 0.0009*(T - 273.15));

    //CuCl-KCl-NaCl with additive
    case 4: return 1000.0 * (3.3112 - 0.0011*(T - 273.15));

    //KNO3–NaNO3–LiNO3–Ca(NO3)2*4H2O 6:1:2:2
    case 5: return 1000.0 * (2.038 - 5.679e-4*(T - 273.15));

    default: throw std::invalid_argument("Invalid molten salt ID");
  }
}

double MoltenSalt::computeConductivity ( double T) {
  switch (_id) {
    
    //Solar salt
    case 0: return 1.16;
    case 1: return 1.16;

    //(LiNaK)2CO3
    case 2: return 0.3362 + 0.000252*(T - 273.15); //Linear regression of data from table 1 of X. An et al.

    //CuCl-KCl-NaCl without additive
    case 3: return 0.4233 + 0.0007*(T - 273.15);
    
    //CuCl-KCl-NaCl with additive
    case 4: return 0.3403 + 0.0008*(T - 273.15);

    //KNO3–NaNO3–LiNO3–Ca(NO3)2*4H2O 6:1:2:2
    case 5: return 0.53;

    default: throw std::invalid_argument("Invalid molten salt ID");
  }
}
