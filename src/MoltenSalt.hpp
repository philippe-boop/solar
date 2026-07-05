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
#ifndef __MOLTEN_SALT_H__
#define __MOLTEN_SALT_H__

#include "Constants.hpp"
#include <cmath>

//Simple class object containing basic thermodynamical properties of molten salt

class MoltenSalt {

private:

  //P.B. 2026-06 : Attribute _viscosity is not necessary since it is only computed at specific temperatures when needed. Attribute _enthalpy is not used either.
  double _temperature; // in K
  double _pressure;    // in kPa
  double _massFlow;    // in kg/s

  int    _id;
  double _meltingPoint;     // in K
  double _maxOperatingTemp; // in K
  double _costPerKg;        // in USD/kg
  double _heatCapacity;     // in J/(kg*K)

public:
  
  MoltenSalt ( double, double, double, int ); //P.B. 2026-06 : Modified this contructor 

  double get_temperature ( void ) const{ return _temperature; }
  double get_pressure    ( void ) const{ return _pressure; }
  double get_massFlow    ( void ) const{ return _massFlow; }

  //P.B. 2026-06
  double get_meltingPoint     ( void ) const{ return _meltingPoint; }
  double get_maxOperatingTemp ( void ) const{ return _maxOperatingTemp;}
  double get_cost             ( void ) const{ return _costPerKg;}
  double get_heatCapacity     ( void ) const{ return _heatCapacity;}

  void set_pressure      ( double pres) { _pressure = pres; }
  void set_massFlow      ( double masf) { _massFlow = masf; }
  void set_temperature   ( double T   ) {_temperature = T;}

  double computeViscosity    ( double ); // Pa*s
  double computeHeatCapacity ( double ); // J/(kg*K)
  double computeDensity      ( double ); // kg/m^3
  double computeConductivity ( double ); // W/(m*K)
};

#endif
