/*
  File:         LogNormal_Distribution.cpp
  Version:      0.0.1
  Date:         03-Feb-2019
  Revision:     03-Feb-2019
  Author:       Jerome Drouin (jerome.p.drouin@gmail.com)

  Editions:	Please go to LogNormal_Distribution.h for Edition Notes.

  LogNormal_Distribution.cpp - Library for 'duino
  https://github.com/newEndeavour/LogNormal_Distribution

  LogNormal_Distribution implements a LogNormal distribution. 

  Copyright (c) 2018-2019 Jerome Drouin  All rights reserved.

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
 
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/


#include "Arduino.h"
#include <math.h>              // required for erf() and M_SQRT1_2
#include "LogNormal_Distribution.h"


// Constructor /////////////////////////////////////////////////////////////////
// Function that handles the creation and setup of instances

LogNormal_Distribution::LogNormal_Distribution(double _Mu, double _Sigma)
{

	// Object parameter's error handling
	error = 1;

	//Std deviation must be positive Real
	if(_Sigma<=0) 	error	= -2; 

	//Set initial values	
	Mu			= _Mu;			// 
	Sigma			= _Sigma;		//

}


// Public Methods //////////////////////////////////////////////////////////////
//Probability Density Function
double LogNormal_Distribution::GetPDF(double x)
{

	//Error check
	if (error<0) return error;
	
	return exp( -0.5 * pow((log(x) - Mu)/Sigma,2) ) / (x * Sigma * sqrt(2.0 * CONSTANT_Pi));
	
}


//Cumulative Distribution Function
double LogNormal_Distribution::GetCDF(double x)
{
double t;
double y;
double z;
double xx = (log(x) - Mu) / Sigma;

	//Error check
	if (error<0) return error;

	return  0.5 * ( 1.0 + erf( M_SQRT1_2 * xx ) );

}


//Mean
double 	LogNormal_Distribution::GetMean(void)
{
	if (error<0)
		return error;
	
	return exp(Mu + 0.5*pow(Sigma,2));

}


//Variance
double 	LogNormal_Distribution::GetVariance(void)
{
	if (error<0)
		return error;

	return exp(pow(Sigma,2)-1) * exp(2.0*Mu + pow(Sigma,2));
}


//Std Deviation
double 	LogNormal_Distribution::GetStdDeviation(void)
{
double var;
double Sigma_sq;
	if (error<0)
		return error;

	Sigma_sq = pow(Sigma,2);
	var = (exp(Sigma_sq)-1) * exp(2.0*Mu + Sigma_sq);
	return sqrt(var);
}


//Skewness
double 	LogNormal_Distribution::GetSkewness(void)
{
double Sigma_sq;

	if (error<0)
		return error;

	Sigma_sq = pow(Sigma,2);
	return (exp(Sigma_sq)+2) * sqrt(exp(Sigma_sq)-1);

}


//Kurtosis
double 	LogNormal_Distribution::GetKurtosis(void)
{
double Sigma_sq;

	if (error<0)
		return error;

	Sigma_sq = pow(Sigma,2); 
	return exp(4.0*Sigma_sq) + 2.0*exp(3.0*Sigma_sq) + 3.0*exp(2.0*Sigma_sq) - 6.0;
}


//Entropy
// log2() Natural Logarithm base 2
double 	LogNormal_Distribution::GetEntropy(void)
{
	if (error<0)
		return error;

	return log2(Sigma * exp(Mu+.5) * sqrt(2.0*CONSTANT_Pi));
}



//Return Quantile z(P) from probability P
/*
			
		        * 
		      *   *
		     *     *
		    *       *
		   *         *
		  *	      *
		*       P       *
       	  |   *  X e [-u.s;+u.s]  *   |
          | *                       * | 
* * * *   |			      |  * * * 
----------+-------------+-------------+--------------------
	-u.s		m	     +u.s


z(P) = LogNormal_Distribution.GetQuantile(P);
ID  P       z(P)
0   0.68    1.000000000000
1   0.80    1.281551565545
2   0.90    1.644853626951  
3   0.95    1.96963984540
4   0.96    2.053749084..
5   0.98    2.326347874041
6   0.99    2.575829303549
7   0.995   2.807033768344
8   0.998   3.090232306168
9   0.999   3.290526731492
10  0.9999  3.890591886413


		        * 
		      *   *
		     *     *
		    *       *
		   *         *
		  *	      *
		*       P       *
       	      *  Xt e [-oo;+t.s]  *   |
            *                       * | 
* * * *    			      |  * * * 
------------------------+-------------+--------------------
-oo			m	     +t.s

z(P) = LogNormal_Distribution.GetQuantile(P);
ID  P       z(P)
0   0.84    1.000000000000
1   0.90    1.281551565545
2   0.95    1.644853626951  
3   0.975   1.96963984540
4   0.980   2.053749084..
5   0.99    2.326347874041
6   0.995   2.575829303549
7   0.9975  2.807033768344
8   0.9990  3.090232306168
9   0.9995  3.290526731492
10  0.99995 3.890591886413

*/
double LogNormal_Distribution::GetQuantile(double p)
{
double Vm;
double Vh = 16;
double Vl = -16;
double Pr;
int i = 0;
double Eps;

	if (error<0)
		return error;

	if (p <= 0.0) {
		return Vl;
	} else if (p >= 1.0) {
        	return Vh;
	} else {        
        	do 
		{
          		i++;
          		Vm = (Vh + Vl) / 2;
            
			Pr = GetCDF(Vm);
          		Eps = abs(Pr - p);
			
			/*
			//DEBUG
			Serial.print("\nF(x) | ");
			Serial.print("i:");
			Serial.print(i);
			Serial.print("\tVl:");
			Serial.print(Vl,4);
			Serial.print("\tVh:");
			Serial.print(Vh,4);
			Serial.print("\tVm:");
			Serial.print(Vm,4);
			Serial.print("\t ->F(Vm):");
			Serial.print(Pr,4);
			Serial.print("\t ->Eps:");
			Serial.print(Eps,4);
			*/
        
          		//New Boundary selection
          		if (Pr > p) {
				Vh = Vm;
				//Serial.print("\t (Pr > p) Vh=Vm->Vh:");
				//Serial.print(Vh,4);
          		} else {
				Vl = Vm;
				//Serial.print("\t (Pr < p) Vl=Vm->Vl:");
				//Serial.print(Vl,4);
			}
            
        	} 
		while ((Eps > CONSTANT_EpsStop) && (i < 70));
	}
            
        if (i >= 70) {
            return -9999;
        } else {
            return Vm;
    	}

}



double LogNormal_Distribution::GetMu(void)
{
	return Mu;
}


double LogNormal_Distribution::GetSigma(void)
{
	return Sigma;
}

void LogNormal_Distribution::SetMu(double _Mu)
{
	//Set initial values	
	Mu			= _Mu;			// 
}


void LogNormal_Distribution::SetSigma(double _Sigma)
{
	// Object parameter's error handling
	error = 1;

	//Std deviation must be positive Real
	if(_Sigma<=0) 	error	= -2; 

	//Set initial values	
	Sigma			= _Sigma;		//
}


// Private Methods /////////////////////////////////////////////////////////////
// Functions only available to other functions in this library


// /////////////////////////////////////////////////////////////////////////////

