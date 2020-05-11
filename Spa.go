package spa

import (
	"errors"
	"math"
	"time"
)

/////////////////////////////////////////////
//                                         //
//      Solar Position Algorithm (SPA)     //
//                   for                   //
//        Solar Radiation Application      //
//                                         //
//               May 12, 2003              //
//                                         //
//                                         //
//   Afshin Michael Andreas                //
//   afshin_andreas@nrel.gov (303)384-6383 //
//                                         //
//   Measurement & Instrumentation Team    //
//   Solar Radiation Research Laboratory   //
//   National Renewable Energy Laboratory  //
//   1617 Cole Blvd, Golden, CO 80401      //
//   This code is based on the NREL        //
//   technical report "Solar Position      //
//   Algorithm for Solar Radiation         //
//   Application" by I. Reda & A. Andreas  //
/////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
//
//   NOTICE
//   Copyright (C) 2008-2011 Alliance for Sustainable Energy, LLC, All Rights Reserved
//
//The Solar Position Algorithm ("Software") is code in development prepared by employees of the
//Alliance for Sustainable Energy, LLC, (hereinafter the "Contractor"), under Contract No.
//DE-AC36-08GO28308 ("Contract") with the U.S. Department of Energy (the "DOE"). The United
//States Government has been granted for itself and others acting on its behalf a paid-up, non-
//exclusive, irrevocable, worldwide license in the Software to reproduce, prepare derivative
//works, and perform publicly and display publicly. Beginning five (5) years after the date
//permission to assert copyright is obtained from the DOE, and subject to any subsequent five
//(5) year renewals, the United States Government is granted for itself and others acting on
//its behalf a paid-up, non-exclusive, irrevocable, worldwide license in the Software to
//reproduce, prepare derivative works, distribute copies to the public, perform publicly and
//display publicly, and to permit others to do so. If the Contractor ceases to make this
//computer software available, it may be obtained from DOE's Office of Scientific and Technical
//Information's Energy Science and Technology Software Center (ESTSC) at P.O. Box 1020, Oak
//Ridge, TN 37831-1020. THIS SOFTWARE IS PROVIDED BY THE CONTRACTOR "AS IS" AND ANY EXPRESS OR
//IMPLIED WARRANTIES, INCLUDING BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
//AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE CONTRACTOR OR THE
//U.S. GOVERNMENT BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
//WHATSOEVER, INCLUDING BUT NOT LIMITED TO CLAIMS ASSOCIATED WITH THE LOSS OF DATA OR PROFITS,
//WHICH MAY RESULT FROM AN ACTION IN CONTRACT, NEGLIGENCE OR OTHER TORTIOUS CLAIM THAT ARISES
//OUT OF OR IN CONNECTION WITH THE ACCESS, USE OR PERFORMANCE OF THIS SOFTWARE.
//
//The Software is being provided for internal, noncommercial purposes only and shall not be
//re-distributed. Please contact Jennifer Ramsey (Jennifer.Ramsey@nrel.gov) in the NREL
//Commercialization and Technology Transfer Office for information concerning a commercial
//license to use the Software.
//
//As a condition of using the Software in an application, the developer of the application
//agrees to reference the use of the Software and make this Notice readily accessible to any
//end-user in a Help|About screen or equivalent manner.
//
///////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////
// Revised 27-FEB-2004 Andreas
//         Added bounds check on inputs and return value for spa_calculate().
// Revised 10-MAY-2004 Andreas
//         Changed temperature bound check minimum from -273.15 to -273 degrees C.
// Revised 17-JUN-2004 Andreas
//         Corrected a problem that caused a bogus sunrise/set/transit on the equinox.
// Revised 18-JUN-2004 Andreas
//         Added a "function" input variable that allows the selecting of desired outputs.
// Revised 21-JUN-2004 Andreas
//         Added 3 new intermediate output values to SPA structure (srha, ssha, & sta).
// Revised 23-JUN-2004 Andreas
//         Enumerations for "function" were renamed and 2 were added.
//         Prevented bound checks on inputs that are not used (based on function).
// Revised 01-SEP-2004 Andreas
//         Changed a local variable from integer to double.
// Revised 12-JUL-2005 Andreas
//         Put a limit on the EOT calculation, so that the result is between -20 and 20.
// Revised 26-OCT-2005 Andreas
//         Set the atmos. refraction correction to zero, when sun is below horizon.
//         Made atmos_refract input a requirement for all "functions".
//         Changed atmos_refract bound check from +/- 10 to +/- 5 degrees.
// Revised 07-NOV-2006 Andreas
//         Corrected 3 earth periodic terms in the L_TERMS array.
//         Corrected 2 earth periodic terms in the R_TERMS array.
// Revised 10-NOV-2006 Andreas
//         Corrected a constant used to calculate topocentric sun declination.
//         Put a limit on observer hour angle, so result is between 0 and 360.
// Revised 13-NOV-2006 Andreas
//         Corrected calculation of topocentric sun declination.
//         Converted all floating point inputs in spa structure to doubles.
// Revised 27-FEB-2007 Andreas
//         Minor correction made as to when atmos. refraction correction is set to zero.
// Revised 21-JAN-2008 Andreas
//         Minor change to two variable declarations.
// Revised 12-JAN-2009 Andreas
//         Changed timezone bound check from +/-12 to +/-18 hours.
// Revised 14-JAN-2009 Andreas
//         Corrected a constant used to calculate ecliptic mean obliquity.
// Revised 01-APR-2013 Andreas
//		   Replace floor with new integer function for tech. report consistency, no affect on results.
//         Add "utility" function prototypes to header file for use with NREL's SAMPA.
//         Rename 4 "utility" function names (remove "sun") for clarity with NREL's SAMPA.
//		   Added delta_ut1 as required input, which the fractional second difference between UT and UTC.
//         Time must be input w/o delta_ut1 adjustment, instead of assuming adjustment was pre-applied.
// Revised 10-JUL-2014 Andreas
//         Change second in spa_data structure from an integer to double to allow fractional second
// Revised 08-SEP-2014 Andreas
//         Corrected description of azm_rotation in header file
//         Limited azimuth180 to range of 0 to 360 deg (instead of -180 to 180) for tech report consistency
//         Changed all variables names from azimuth180 to azimuth_astro
//         Renamed 2 "utility" function names for consistency
///////////////////////////////////////////////////////////////////////////////////////////////

const SunRadius float64 = 0.26667

const LCount int64 = 6
const BCount int64 = 2
const RCount int64 = 5
const YCount int64 = 63

const LMaxSubcount int64 = 64
const BMaxSubcount int64 = 5
const RMaxSubcount int64 = 40

const TermYCount int64 = TermXCount

var lSubcount = []int64{64, 34, 20, 7, 3, 1}
var bSubcount = []int64{5, 2}
var rSubcount = []int64{40, 10, 6, 2, 1}

const (
	TermA     = 0
	TermB     = 1
	TermC     = 2
	TermCount = 3
)
const (
	TermX0     = 0
	TermX1     = 1
	TermX2     = 2
	TermX3     = 3
	TermX4     = 4
	TermXCount = 5
)

const (
	TermPsiA    = 0
	TermPsiB    = 1
	TermEpsC    = 2
	TermEpsD    = 3
	TermPeCount = 4
)
const (
	JdMinus = 0
	JdZero  = 1
	JdPlus  = 2
	JdCount = 3
)
const (
	SunTransit = 0
	SunRise    = 1
	SunSet     = 2
	SunCount   = 3
)

///////////////////////////////////////////////////
///  Earth Periodic Terms
///////////////////////////////////////////////////
var LTerms = [][][]float64{{
	{175347046.0, 0, 0},
	{3341656.0, 4.6692568, 6283.07585},
	{34894.0, 4.6261, 12566.1517},
	{3497.0, 2.7441, 5753.3849},
	{3418.0, 2.8289, 3.5231},
	{3136.0, 3.6277, 77713.7715},
	{2676.0, 4.4181, 7860.4194},
	{2343.0, 6.1352, 3930.2097},
	{1324.0, 0.7425, 11506.7698},
	{1273.0, 2.0371, 529.691},
	{1199.0, 1.1096, 1577.3435},
	{990, 5.233, 5884.927},
	{902, 2.045, 26.298},
	{857, 3.508, 398.149},
	{780, 1.179, 5223.694},
	{753, 2.533, 5507.553},
	{505, 4.583, 18849.228},
	{492, 4.205, 775.523},
	{357, 2.92, 0.067},
	{317, 5.849, 11790.629},
	{284, 1.899, 796.298},
	{271, 0.315, 10977.079},
	{243, 0.345, 5486.778},
	{206, 4.806, 2544.314},
	{205, 1.869, 5573.143},
	{202, 2.458, 6069.777},
	{156, 0.833, 213.299},
	{132, 3.411, 2942.463},
	{126, 1.083, 20.775},
	{115, 0.645, 0.98},
	{103, 0.636, 4694.003},
	{102, 0.976, 15720.839},
	{102, 4.267, 7.114},
	{99, 6.21, 2146.17},
	{98, 0.68, 155.42},
	{86, 5.98, 161000.69},
	{85, 1.3, 6275.96},
	{85, 3.67, 71430.7},
	{80, 1.81, 17260.15},
	{79, 3.04, 12036.46},
	{75, 1.76, 5088.63},
	{74, 3.5, 3154.69},
	{74, 4.68, 801.82},
	{70, 0.83, 9437.76},
	{62, 3.98, 8827.39},
	{61, 1.82, 7084.9},
	{57, 2.78, 6286.6},
	{56, 4.39, 14143.5},
	{56, 3.47, 6279.55},
	{52, 0.19, 12139.55},
	{52, 1.33, 1748.02},
	{51, 0.28, 5856.48},
	{49, 0.49, 1194.45},
	{41, 5.37, 8429.24},
	{41, 2.4, 19651.05},
	{39, 6.17, 10447.39},
	{37, 6.04, 10213.29},
	{37, 2.57, 1059.38},
	{36, 1.71, 2352.87},
	{36, 1.78, 6812.77},
	{33, 0.59, 17789.85},
	{30, 0.44, 83996.85},
	{30, 2.74, 1349.87},
	{25, 3.16, 4690.48}}, {
	{628331966747.0, 0, 0},
	{206059.0, 2.678235, 6283.07585},
	{4303.0, 2.6351, 12566.1517},
	{425.0, 1.59, 3.523},
	{119.0, 5.796, 26.298},
	{109.0, 2.966, 1577.344},
	{93, 2.59, 18849.23},
	{72, 1.14, 529.69},
	{68, 1.87, 398.15},
	{67, 4.41, 5507.55},
	{59, 2.89, 5223.69},
	{56, 2.17, 155.42},
	{45, 0.4, 796.3},
	{36, 0.47, 775.52},
	{29, 2.65, 7.11},
	{21, 5.34, 0.98},
	{19, 1.85, 5486.78},
	{19, 4.97, 213.3},
	{17, 2.99, 6275.96},
	{16, 0.03, 2544.31},
	{16, 1.43, 2146.17},
	{15, 1.21, 10977.08},
	{12, 2.83, 1748.02},
	{12, 3.26, 5088.63},
	{12, 5.27, 1194.45},
	{12, 2.08, 4694},
	{11, 0.77, 553.57},
	{10, 1.3, 6286.6},
	{10, 4.24, 1349.87},
	{9, 2.7, 242.73},
	{9, 5.64, 951.72},
	{8, 5.3, 2352.87},
	{6, 2.65, 9437.76},
	{6, 4.67, 4690.48}}, {
	{52919.0, 0, 0},
	{8720.0, 1.0721, 6283.0758},
	{309.0, 0.867, 12566.152},
	{27, 0.05, 3.52},
	{16, 5.19, 26.3},
	{16, 3.68, 155.42},
	{10, 0.76, 18849.23},
	{9, 2.06, 77713.77},
	{7, 0.83, 775.52},
	{5, 4.66, 1577.34},
	{4, 1.03, 7.11},
	{4, 3.44, 5573.14},
	{3, 5.14, 796.3},
	{3, 6.05, 5507.55},
	{3, 1.19, 242.73},
	{3, 6.12, 529.69},
	{3, 0.31, 398.15},
	{3, 2.28, 553.57},
	{2, 4.38, 5223.69},
	{2, 3.75, 0.98}},
	{
		{289.0, 5.844, 6283.076},
		{35, 0, 0},
		{17, 5.49, 12566.15},
		{3, 5.2, 155.42},
		{1, 4.72, 3.52},
		{1, 5.3, 18849.23},
		{1, 5.97, 242.73}},
	{
		{114.0, 3.142, 0},
		{8, 4.13, 6283.08},
		{1, 3.84, 12566.15}},
	{{1, 3.14, 0}}}

var BTerms = [][][]float64{{
	{280.0, 3.199, 84334.662},
	{102.0, 5.422, 5507.553},
	{80, 3.88, 5223.69},
	{44, 3.7, 2352.87},
	{32, 4, 1577.34}}, {
	{9, 3.9, 5507.55},
	{6, 1.73, 5223.69}}}

var RTerms = [][][]float64{{{100013989.0, 0, 0},
	{1670700.0, 3.0984635, 6283.07585},
	{13956.0, 3.05525, 12566.1517},
	{3084.0, 5.1985, 77713.7715},
	{1628.0, 1.1739, 5753.3849},
	{1576.0, 2.8469, 7860.4194},
	{925.0, 5.453, 11506.77},
	{542.0, 4.564, 3930.21},
	{472.0, 3.661, 5884.927},
	{346.0, 0.964, 5507.553},
	{329.0, 5.9, 5223.694},
	{307.0, 0.299, 5573.143},
	{243.0, 4.273, 11790.629},
	{212.0, 5.847, 1577.344},
	{186.0, 5.022, 10977.079},
	{175.0, 3.012, 18849.228},
	{110.0, 5.055, 5486.778},
	{98, 0.89, 6069.78},
	{86, 5.69, 15720.84},
	{86, 1.27, 161000.69},
	{65, 0.27, 17260.15},
	{63, 0.92, 529.69},
	{57, 2.01, 83996.85},
	{56, 5.24, 71430.7},
	{49, 3.25, 2544.31},
	{47, 2.58, 775.52},
	{45, 5.54, 9437.76},
	{43, 6.01, 6275.96},
	{39, 5.36, 4694},
	{38, 2.39, 8827.39},
	{37, 0.83, 19651.05},
	{37, 4.9, 12139.55},
	{36, 1.67, 12036.46},
	{35, 1.84, 2942.46},
	{33, 0.24, 7084.9},
	{32, 0.18, 5088.63},
	{32, 1.78, 398.15},
	{28, 1.21, 6286.6},
	{28, 1.9, 6279.55},
	{26, 4.59, 10447.39}},
	{
		{103019.0, 1.10749, 6283.07585},
		{1721.0, 1.0644, 12566.1517},
		{702.0, 3.142, 0},
		{32, 1.02, 18849.23},
		{31, 2.84, 5507.55},
		{25, 1.32, 5223.69},
		{18, 1.42, 1577.34},
		{10, 5.91, 10977.08},
		{9, 1.42, 6275.96},
		{9, 0.27, 5486.78}},
	{
		{4359.0, 5.7846, 6283.0758},
		{124.0, 5.579, 12566.152},
		{12, 3.14, 0},
		{9, 3.63, 77713.77},
		{6, 1.87, 5573.14},
		{3, 5.47, 18849.23}},
	{
		{145.0, 4.273, 6283.076},
		{7, 3.92, 12566.15}},
	{
		{4, 2.56, 6283.08}}}

////////////////////////////////////////////////////////////////
///  Periodic Terms for the nutation in longitude and obliquity
////////////////////////////////////////////////////////////////

var YTerms = [][]int64{
	{0, 0, 0, 0, 1},
	{-2, 0, 0, 2, 2},
	{0, 0, 0, 2, 2},
	{0, 0, 0, 0, 2},
	{0, 1, 0, 0, 0},
	{0, 0, 1, 0, 0},
	{-2, 1, 0, 2, 2},
	{0, 0, 0, 2, 1},
	{0, 0, 1, 2, 2},
	{-2, -1, 0, 2, 2},
	{-2, 0, 1, 0, 0},
	{-2, 0, 0, 2, 1},
	{0, 0, -1, 2, 2},
	{2, 0, 0, 0, 0},
	{0, 0, 1, 0, 1},
	{2, 0, -1, 2, 2},
	{0, 0, -1, 0, 1},
	{0, 0, 1, 2, 1},
	{-2, 0, 2, 0, 0},
	{0, 0, -2, 2, 1},
	{2, 0, 0, 2, 2},
	{0, 0, 2, 2, 2},
	{0, 0, 2, 0, 0},
	{-2, 0, 1, 2, 2},
	{0, 0, 0, 2, 0},
	{-2, 0, 0, 2, 0},
	{0, 0, -1, 2, 1},
	{0, 2, 0, 0, 0},
	{2, 0, -1, 0, 1},
	{-2, 2, 0, 2, 2},
	{0, 1, 0, 0, 1},
	{-2, 0, 1, 0, 1},
	{0, -1, 0, 0, 1},
	{0, 0, 2, -2, 0},
	{2, 0, -1, 2, 1},
	{2, 0, 1, 2, 2},
	{0, 1, 0, 2, 2},
	{-2, 1, 1, 0, 0},
	{0, -1, 0, 2, 2},
	{2, 0, 0, 2, 1},
	{2, 0, 1, 0, 0},
	{-2, 0, 2, 2, 2},
	{-2, 0, 1, 2, 1},
	{2, 0, -2, 0, 1},
	{2, 0, 0, 0, 1},
	{0, -1, 1, 0, 0},
	{-2, -1, 0, 2, 1},
	{-2, 0, 0, 0, 1},
	{0, 0, 2, 2, 1},
	{-2, 0, 2, 0, 1},
	{-2, 1, 0, 2, 1},
	{0, 0, 1, -2, 0},
	{-1, 0, 1, 0, 0},
	{-2, 1, 0, 0, 0},
	{1, 0, 0, 0, 0},
	{0, 0, 1, 2, 0},
	{0, 0, -2, 2, 2},
	{-1, -1, 1, 0, 0},
	{0, 1, 1, 0, 0},
	{0, -1, 1, 2, 2},
	{2, -1, -1, 2, 2},
	{0, 0, 3, 2, 2},
	{2, -1, 0, 2, 2},
}
var PeTerms = [][]float64{
	{-171996, -174.2, 92025, 8.9},
	{-13187, -1.6, 5736, -3.1},
	{-2274, -0.2, 977, -0.5},
	{2062, 0.2, -895, 0.5},
	{1426, -3.4, 54, -0.1},
	{712, 0.1, -7, 0},
	{-517, 1.2, 224, -0.6},
	{-386, -0.4, 200, 0},
	{-301, 0, 129, -0.1},
	{217, -0.5, -95, 0.3},
	{-158, 0, 0, 0},
	{129, 0.1, -70, 0},
	{123, 0, -53, 0},
	{63, 0, 0, 0},
	{63, 0.1, -33, 0},
	{-59, 0, 26, 0},
	{-58, -0.1, 32, 0},
	{-51, 0, 27, 0},
	{48, 0, 0, 0},
	{46, 0, -24, 0},
	{-38, 0, 16, 0},
	{-31, 0, 13, 0},
	{29, 0, 0, 0},
	{29, 0, -12, 0},
	{26, 0, 0, 0},
	{-22, 0, 0, 0},
	{21, 0, -10, 0},
	{17, -0.1, 0, 0},
	{16, 0, -8, 0},
	{-16, 0.1, 7, 0},
	{-15, 0, 9, 0},
	{-13, 0, 7, 0},
	{-12, 0, 6, 0},
	{11, 0, 0, 0},
	{-10, 0, 5, 0},
	{-8, 0, 3, 0},
	{7, 0, -3, 0},
	{-7, 0, 0, 0},
	{-7, 0, 3, 0},
	{-7, 0, 3, 0},
	{6, 0, 0, 0},
	{6, 0, -3, 0},
	{6, 0, -3, 0},
	{-6, 0, 3, 0},
	{-6, 0, 3, 0},
	{5, 0, 0, 0},
	{-5, 0, 3, 0},
	{-5, 0, 3, 0},
	{-5, 0, 3, 0},
	{4, 0, 0, 0},
	{4, 0, 0, 0},
	{4, 0, 0, 0},
	{-4, 0, 0, 0},
	{-4, 0, 0, 0},
	{-4, 0, 0, 0},
	{3, 0, 0, 0},
	{-3, 0, 0, 0},
	{-3, 0, 0, 0},
	{-3, 0, 0, 0},
	{-3, 0, 0, 0},
	{-3, 0, 0, 0},
	{-3, 0, 0, 0},
	{-3, 0, 0, 0},
}

///////////////////////////////////////////////

// Spa interface defines the public functions
type Spa interface {
	Calculate() error
	//-----------------INPUTE VALUES--------------------
	// Helper function to use date
	SetDate(time time.Time)
	GetDate() time.Time
	// 4-digit year,      valid range: -2000 to 6000
	SetYear(int)
	GetYear() int
	// 2-digit month,         valid range: 1 to  12
	SetMonth(int)
	GetMonth() int
	// Observer local hour,   valid range: 0 to  24
	SetDay(int)
	GetDay() int
	// Observer local minute, valid range: 0 to  59
	SetMinute(int)
	GetMinute() int
	// Observer local second, valid range: 0 to <60
	SetSecond(float64)
	GetSecond() float64
	// Fractional second difference between UTC and UT which is used to adjust UTC for earth's irregular rotation rate and is derived
	// from observation only and is reported in this bulletin: http://maia.usno.navy.mil/ser7/ser7.dat where delta_ut1 = DUT1
	// valid range: -1 to 1 second (exclusive)
	SetDeltaUt1(float64)
	GetDeltaUt1() float64

	// Difference between earth rotation time and terrestrial time It is derived from observation only and is reported in this
	// bulletin: http://maia.usno.navy.mil/ser7/ser7.dat, where delta_t = 32.184 + (TAI-UTC) - DUT1
	// valid range: -8000 to 8000 seconds
	SetDeltaT(float64)
	GetDeltaT() float64
	// Observer time zone (negative west of Greenwich) valid range: -18   to   18 hours
	SetTimezone(float64)
	GetTimezone() float64
	// Observer longitude (negative west of Greenwich) valid range: -180  to  180 degrees
	SetLongitude(float64)
	GetLongitude() float64
	// Observer latitude (negative south of equator) valid range: -90   to   90 degrees
	SetLatitude(float64)
	GetLatitude() float64
	// Observer elevation [meters] valid range: -6500000 or higher meters
	SetElevation(float64)
	GetElevation() float64
	// Annual average local pressure [millibars] valid range:    0 to 5000 millibars
	SetPressure(float64)
	GetPressure() float64
	// Annual average local temperature [degrees Celsius] valid range: -273 to 6000 degrees Celsius
	SetTemperature(float64)
	GetTemperature() float64
	// Surface slope (measured from the horizontal plane) valid range: -360 to 360 degrees
	SetSlope(float64)
	GetSlope() float64
	// Surface azimuth rotation (measured from south to projection of surface normal on horizontal plane, negative east) valid range: -360 to 360 degrees
	SetAzmRotation(float64)
	GetAzmRotation() float64
	// Atmospheric refraction at sunrise and sunset (0.5667 deg is typical)valid range: -5   to   5 degrees
	SetAtmosRefract(float64)
	GetAtmosRefract() float64
	// Switch to choose functions for desired output (from enumeration)
	SetSPAFunction(SPAFunctions)
	GetSPAFunction() SPAFunctions
	//-----------------Intermediate OUTPUT VALUES--------------------
	//Julian day
	GetJd() float64
	//Julian century
	GetJc() float64
	//Julian ephemeris day
	GetJde() float64
	//Julian ephemeris century
	GetJce() float64
	//Julian ephemeris millennium
	GetJme() float64
	//earth heliocentric longitude [degrees]
	GetL() float64
	//earth heliocentric latitude [degrees]
	GetB() float64
	//earth radius vector [Astronomical Units, AU]
	GetR() float64
	//geocentric longitude [degrees]
	GetTheta() float64
	//geocentric latitude [degrees]
	GetBeta() float64
	//mean elongation (moon-sun) [degrees]
	GetX0() float64
	//mean anomaly (sun) [degrees]
	GetX1() float64
	//mean anomaly (moon) [degrees]
	GetX2() float64
	//argument latitude (moon) [degrees]
	GetX3() float64
	//ascending longitude (moon) [degrees]
	GetX4() float64
	//nutation longitude [degrees]
	GetDelPsi() float64
	//nutation obliquity [degrees]
	GetDelEpsilon() float64
	//ecliptic mean obliquity [arc seconds]
	GetEpsilon0() float64
	//ecliptic true obliquity  [degrees]
	GetEpsilon() float64
	//aberration correction [degrees]
	GetDelTau() float64
	//apparent sun longitude [degrees]
	GetLamda() float64
	//Greenwich mean sidereal time [degrees]
	GetNu0() float64
	//Greenwich sidereal time [degrees]
	GetNu() float64
	//geocentric sun right ascension [degrees]
	GetAlpha() float64
	//geocentric sun declination [degrees]
	GetDelta() float64
	//observer hour angle [degrees]
	GetH() float64
	//sun equatorial horizontal parallax [degrees]
	GetXi() float64
	//sun right ascension parallax [degrees]
	GetDelAlpha() float64
	//topocentric sun declination [degrees]
	GetDeltaPrime() float64
	//topocentric sun right ascension [degrees]
	GetAlphaPrime() float64
	//topocentric local hour angle [degrees]
	GetHPrime() float64
	//topocentric elevation angle (uncorrected) [degrees]
	GetE0() float64
	//atmospheric refraction correction [degrees]
	GetDelE() float64
	//topocentric elevation angle (corrected) [degrees]
	GetE() float64
	//equation of time [minutes]
	GetEot() float64
	//sunrise hour angle [degrees]
	GetSrha() float64
	//sunset hour angle [degrees]
	GetSsha() float64
	//sun transit altitude [degrees]
	GetSta() float64
	//---------------------Final OUTPUT VALUES------------------------
	//topocentric zenith angle [degrees]
	GetZenith() float64
	//topocentric azimuth angle (westward from south) [for astronomers]
	GetAzimuthAstro() float64
	//topocentric azimuth angle (eastward from north) [for navigators and solar radiation]
	GetAzimuth() float64
	//surface incidence angle [degrees]
	GetIncidence() float64
	//local sun transit time (or solar noon) [fractional hour]
	GetSuntransit() float64
	//local sunrise time (+/- 30 seconds) [fractional hour]
	GetSunrise() time.Time
	//local sunset time (+/- 30 seconds) [fractional hour]
	GetSunset() time.Time
}

// NewSpa creates new SPA instance
func NewSpa(dt time.Time, latitude float64, longitude float64, elevation float64, pressure float64, temperature float64, deltaT float64, deltaUt1 float64, slope float64, azmRotation float64, atmosRefract float64) (Spa, error) {
	var s spa
	s.init()
	s.SetDate(dt)
	s.latitude = latitude
	s.longitude = longitude
	s.elevation = elevation
	s.pressure = pressure
	s.temperature = temperature
	s.deltaT = deltaT
	s.deltaUt1 = deltaUt1
	s.slope = slope
	s.azmRotation = azmRotation
	s.atmosRefract = atmosRefract
	s.function = SpaAll

	return &s, s.Calculate()
}

type spa struct {
	//----------------------INPUT VALUES------------------------

	year   int     // 4-digit year,      valid range: -2000 to 6000, error code: 1
	month  int     // 2-digit month,         valid range: 1 to  12,  error code: 2
	day    int     // 2-digit day,           valid range: 1 to  31,  error code: 3
	hour   int     // Observer local hour,   valid range: 0 to  24,  error code: 4
	minute int     // Observer local minute, valid range: 0 to  59,  error code: 5
	second float64 // Observer local second, valid range: 0 to <60,  error code: 6

	deltaUt1 float64 // Fractional second difference between UTC and UT which is used
	// to adjust UTC for earth's irregular rotation rate and is derived
	// from observation only and is reported in this bulletin:
	// http://maia.usno.navy.mil/ser7/ser7.dat,
	// where delta_ut1 = DUT1
	// valid range: -1 to 1 second (exclusive), error code 17

	deltaT float64 // Difference between earth rotation time and terrestrial time
	// It is derived from observation only and is reported in this
	// bulletin: http://maia.usno.navy.mil/ser7/ser7.dat,
	// where delta_t = 32.184 + (TAI-UTC) - DUT1
	// valid range: -8000 to 8000 seconds, error code: 7

	timezone float64 // Observer time zone (negative west of Greenwich)
	// valid range: -18   to   18 hours,   error code: 8

	longitude float64 // Observer longitude (negative west of Greenwich)
	// valid range: -180  to  180 degrees, error code: 9

	latitude float64 // Observer latitude (negative south of equator)
	// valid range: -90   to   90 degrees, error code: 10

	elevation float64 // Observer elevation [meters]
	// valid range: -6500000 or higher meters,    error code: 11

	pressure float64 // Annual average local pressure [millibars]
	// valid range:    0 to 5000 millibars,       error code: 12

	temperature float64 // Annual average local temperature [degrees Celsius]
	// valid range: -273 to 6000 degrees Celsius, error code; 13

	slope float64 // Surface slope (measured from the horizontal plane)
	// valid range: -360 to 360 degrees, error code: 14

	azmRotation float64 // Surface azimuth rotation (measured from south to projection of
	//     surface normal on horizontal plane, negative east)
	// valid range: -360 to 360 degrees, error code: 15

	atmosRefract float64 // Atmospheric refraction at sunrise and sunset (0.5667 deg is typical)
	// valid range: -5   to   5 degrees, error code: 16

	function SPAFunctions // Switch to choose functions for desired output (from enumeration)

	//-----------------Intermediate OUTPUT VALUES--------------------

	jd float64 //Julian day
	jc float64 //Julian century

	jde float64 //Julian ephemeris day
	jce float64 //Julian ephemeris century
	jme float64 //Julian ephemeris millennium

	l float64 //earth heliocentric longitude [degrees]
	b float64 //earth heliocentric latitude [degrees]
	r float64 //earth radius vector [Astronomical Units, AU]

	theta float64 //geocentric longitude [degrees]
	beta  float64 //geocentric latitude [degrees]

	x0 float64 //mean elongation (moon-sun) [degrees]
	x1 float64 //mean anomaly (sun) [degrees]
	x2 float64 //mean anomaly (moon) [degrees]
	x3 float64 //argument latitude (moon) [degrees]
	x4 float64 //ascending longitude (moon) [degrees]

	delPsi     float64 //nutation longitude [degrees]
	delEpsilon float64 //nutation obliquity [degrees]
	epsilon0   float64 //ecliptic mean obliquity [arc seconds]
	epsilon    float64 //ecliptic true obliquity  [degrees]

	delTau float64 //aberration correction [degrees]
	lamda  float64 //apparent sun longitude [degrees]
	nu0    float64 //Greenwich mean sidereal time [degrees]
	nu     float64 //Greenwich sidereal time [degrees]

	alpha float64 //geocentric sun right ascension [degrees]
	delta float64 //geocentric sun declination [degrees]

	h          float64 //observer hour angle [degrees]
	xi         float64 //sun equatorial horizontal parallax [degrees]
	delAlpha   float64 //sun right ascension parallax [degrees]
	deltaPrime float64 //topocentric sun declination [degrees]
	alphaPrime float64 //topocentric sun right ascension [degrees]
	hPrime     float64 //topocentric local hour angle [degrees]

	e0   float64 //topocentric elevation angle (uncorrected) [degrees]
	delE float64 //atmospheric refraction correction [degrees]
	e    float64 //topocentric elevation angle (corrected) [degrees]

	eot  float64 //equation of time [minutes]
	srha float64 //sunrise hour angle [degrees]
	ssha float64 //sunset hour angle [degrees]
	sta  float64 //sun transit altitude [degrees]
	mRts []float64

	//---------------------Final OUTPUT VALUES------------------------

	zenith       float64 //topocentric zenith angle [degrees]
	azimuthAstro float64 //topocentric azimuth angle (westward from south) [for astronomers]
	azimuth      float64 //topocentric azimuth angle (eastward from north) [for navigators and solar radiation]
	incidence    float64 //surface incidence angle [degrees]

	suntransit float64 //local sun transit time (or solar noon) [fractional hour]
	sunrise    float64 //local sunrise time (+/- 30 seconds) [fractional hour]
	sunset     float64 //local sunset time (+/- 30 seconds) [fractional hour]

}

func (s *spa) GetJd() float64 {
	return s.jd
}

func (s *spa) GetJc() float64 {
	return s.jc
}

func (s *spa) GetJde() float64 {
	return s.jde
}

func (s *spa) GetJce() float64 {
	return s.jce
}

func (s *spa) GetJme() float64 {
	return s.jme
}

func (s *spa) GetL() float64 {
	return s.l
}

func (s *spa) GetB() float64 {
	return s.b
}

func (s *spa) GetR() float64 {
	return s.r
}

func (s *spa) GetTheta() float64 {
	return s.theta
}

func (s *spa) GetBeta() float64 {
	return s.beta
}

func (s *spa) GetX0() float64 {
	return s.x0
}

func (s *spa) GetX1() float64 {
	return s.x1
}

func (s *spa) GetX2() float64 {
	return s.x2
}

func (s *spa) GetX3() float64 {
	return s.x3
}

func (s *spa) GetX4() float64 {
	return s.x4
}

func (s *spa) GetDelPsi() float64 {
	return s.delPsi
}

func (s *spa) GetDelEpsilon() float64 {
	return s.delEpsilon
}

func (s *spa) GetEpsilon0() float64 {
	return s.epsilon0
}

func (s *spa) GetEpsilon() float64 {
	return s.epsilon
}

func (s *spa) GetDelTau() float64 {
	return s.delTau
}

func (s *spa) GetLamda() float64 {
	return s.lamda
}

func (s *spa) GetNu0() float64 {
	return s.nu0
}

func (s *spa) GetNu() float64 {
	return s.nu
}

func (s *spa) GetAlpha() float64 {
	return s.alpha
}

func (s *spa) GetDelta() float64 {
	return s.delta
}

func (s *spa) GetH() float64 {
	return s.h
}

func (s *spa) GetXi() float64 {
	return s.xi
}

func (s *spa) GetDelAlpha() float64 {
	return s.delAlpha
}

func (s *spa) GetDeltaPrime() float64 {
	return s.deltaPrime
}

func (s *spa) GetAlphaPrime() float64 {
	return s.alphaPrime
}

func (s *spa) GetHPrime() float64 {
	return s.hPrime
}

func (s *spa) GetE0() float64 {
	return s.e0
}

func (s *spa) GetDelE() float64 {
	return s.delE
}

func (s *spa) GetE() float64 {
	return s.e
}

func (s *spa) GetEot() float64 {
	return s.eot
}

func (s *spa) GetSrha() float64 {
	return s.srha
}

func (s *spa) GetSsha() float64 {
	return s.ssha
}

func (s *spa) GetSta() float64 {
	return s.sta
}

func (s *spa) GetZenith() float64 {
	return s.zenith
}

func (s *spa) GetAzimuthAstro() float64 {
	return s.azimuthAstro
}

func (s *spa) GetAzimuth() float64 {
	return s.azimuth
}

func (s *spa) GetIncidence() float64 {
	return s.incidence
}

func (s *spa) GetSuntransit() float64 {
	return s.suntransit
}

func (s *spa) GetSunrise() time.Time {
	h, m, sec := s.calculateHourMinSec(s.sunrise)
	dt := time.Date(s.year, time.Month(s.month), s.day, 0, 0, 0, 0, time.FixedZone("ManualTimeZone", int(s.timezone*3600)))
	return dt.Add(time.Hour*time.Duration(h) +
		time.Minute*time.Duration(m) +
		time.Second*time.Duration(sec))
}

func (s *spa) GetSunset() time.Time {
	h, m, sec := s.calculateHourMinSec(s.sunset)
	dt := time.Date(s.year, time.Month(s.month), s.day, 0, 0, 0, 0, time.FixedZone("ManualTimeZone", int(s.timezone*3600)))
	return dt.Add(time.Hour*time.Duration(h) +
		time.Minute*time.Duration(m) +
		time.Second*time.Duration(sec))
}

func (s *spa) SetDate(dt time.Time) {
	_, offset := dt.Zone()
	s.year = dt.Year()
	s.month = int(dt.Month())
	s.day = dt.Day()
	s.hour = dt.Hour()
	s.minute = dt.Minute()
	s.second = float64(dt.Second())
	s.timezone = float64(offset / 3600)

}

func (s *spa) GetDate() time.Time {
	return time.Date(s.year, time.Month(s.month), s.day, s.hour, s.minute, int(s.second), 0, time.FixedZone("ManualTimeZone", int(s.timezone*3600)))
}
func (s *spa) calculateHourMinSec(decHours float64) (hours int, minutes int, seconds int) {
	min := 60.0 * (decHours - float64(int(decHours)))
	sec := 60.0 * (min - float64(int(min)))

	return int(math.Floor(decHours)), int(min), int(sec)
}
func (s *spa) SetYear(year int) {
	s.year = year
}

func (s *spa) GetYear() int {
	return s.year
}

func (s *spa) SetMonth(month int) {
	s.month = month
}

func (s *spa) GetMonth() int {
	return s.month
}

func (s *spa) SetDay(day int) {
	s.day = day
}

func (s *spa) GetDay() int {
	return s.day
}

func (s *spa) SetMinute(minute int) {
	s.minute = minute
}

func (s *spa) GetMinute() int {
	return s.minute
}

func (s *spa) SetSecond(second float64) {
	s.second = second
}

func (s *spa) GetSecond() float64 {
	return s.second
}

func (s *spa) SetDeltaUt1(deltaUt1 float64) {
	s.deltaUt1 = deltaUt1
}

func (s *spa) GetDeltaUt1() float64 {
	return s.deltaUt1
}

func (s *spa) SetDeltaT(deltaT float64) {
	s.deltaT = deltaT
}

func (s *spa) GetDeltaT() float64 {
	return s.deltaT
}

func (s *spa) SetTimezone(tz float64) {
	s.timezone = tz
}

func (s *spa) GetTimezone() float64 {
	return s.timezone
}

func (s *spa) SetLongitude(lon float64) {
	s.longitude = lon
}

func (s *spa) GetLongitude() float64 {
	return s.longitude
}

func (s *spa) SetLatitude(lat float64) {
	s.latitude = lat
}

func (s *spa) GetLatitude() float64 {
	return s.latitude
}

func (s *spa) SetElevation(elevation float64) {
	s.elevation = elevation
}

func (s *spa) GetElevation() float64 {
	return s.elevation
}

func (s *spa) SetPressure(pressure float64) {
	s.pressure = pressure
}

func (s *spa) GetPressure() float64 {
	return s.pressure
}

func (s *spa) SetTemperature(temp float64) {
	s.temperature = temp
}

func (s *spa) GetTemperature() float64 {
	return s.temperature
}

func (s *spa) SetSlope(slope float64) {
	s.slope = slope
}

func (s *spa) GetSlope() float64 {
	return s.slope
}

func (s *spa) SetAzmRotation(azmRotation float64) {
	s.azmRotation = azmRotation
}

func (s *spa) GetAzmRotation() float64 {
	return s.azmRotation
}

func (s *spa) SetAtmosRefract(atmosRefract float64) {
	s.atmosRefract = atmosRefract
}

func (s *spa) GetAtmosRefract() float64 {
	return s.atmosRefract
}

func (s *spa) SetSPAFunction(functions SPAFunctions) {
	s.function = functions
}

func (s *spa) GetSPAFunction() SPAFunctions {
	return s.function
}

func (s *spa) init() {
	// use  some dummy values for init
	s.year = 2003
	s.month = 10
	s.day = 17
	s.hour = 12
	s.minute = 30
	s.second = 30
	s.timezone = -7.0
	s.deltaUt1 = 0
	s.deltaT = 67
	s.longitude = -105.1786
	s.latitude = 39.742476
	s.elevation = 1830.14
	s.pressure = 820
	s.temperature = 11
	s.slope = 30
	s.azmRotation = -10
	s.atmosRefract = 0.5667
	s.function = SpaAll
}

//Calculate SPA output values (in structure) based on input values passed in structure
func (s *spa) Calculate() error {

	// renew the date
	s.SetDate(s.GetDate())

	err := s.validate()
	if err != nil {
		return err
	}

	s.jd = s.julianDay(s.year, s.month, s.day, s.hour,
		s.minute, s.second, s.deltaUt1, s.timezone)

	s.calculateGeocentricSunRightAscensionAndDeclination()

	s.h = s.observerHourAngle(s.nu, s.longitude, s.alpha)
	s.xi = s.sunEquatorialHorizontalParallax(s.r)

	s.rightAscensionParallaxAndTopocentricDec(s.latitude, s.elevation, s.xi, s.h, s.delta)

	s.alphaPrime = s.topocentricRightAscension(s.alpha, s.delAlpha)
	s.hPrime = s.topocentricLocalHourAngle(s.h, s.delAlpha)

	s.e0 = s.topocentricElevationAngle(s.latitude, s.deltaPrime, s.hPrime)
	s.delE = s.atmosphericRefractionCorrection(s.pressure, s.temperature,
		s.atmosRefract, s.e0)
	s.e = s.topocentricElevationAngleCorrected(s.e0, s.delE)

	s.zenith = s.topocentricZenithAngle(s.e)
	s.azimuthAstro = s.topocentricAzimuthAngleAstro(s.hPrime, s.latitude,
		s.deltaPrime)
	s.azimuth = s.topocentricAzimuthAngle(s.azimuthAstro)

	if (s.function == SpaZaInc) || (s.function == SpaAll) {
		s.incidence = s.surfaceIncidenceAngle(s.zenith, s.azimuthAstro,
			s.azmRotation, s.slope)
	}
	if (s.function == SpaZaRts) || (s.function == SpaAll) {
		s.calculateEotAndSunRiseTransitSet()
	}

	return nil
}
func (s *spa) deg2rad(degrees float64) float64 {
	return (math.Pi / 180.0) * degrees
}
func (s *spa) rad2deg(radians float64) float64 {
	return (180.0 / math.Pi) * radians
}
func (s *spa) limitDegrees(degrees float64) float64 {
	var limited float64
	degrees /= 360.0
	limited = 360.0 * (degrees - math.Floor(degrees))
	if limited < 0 {
		limited += 360.0
	}

	return limited
}
func (s *spa) limitDegrees180pm(degrees float64) float64 {
	var limited float64

	degrees /= 360.0
	limited = 360.0 * (degrees - math.Floor(degrees))
	if limited < -180.0 {
		limited += 360.0
	} else if limited > 180.0 {
		limited -= 360.0
	}

	return limited
}

func (s *spa) limitDegrees180(degrees float64) float64 {
	var limited float64

	degrees /= 180.0
	limited = 180.0 * (degrees - math.Floor(degrees))
	if limited < 0 {
		limited += 180.0
	}

	return limited
}

func (s *spa) limitZero2one(value float64) float64 {
	var limited float64

	limited = value - math.Floor(value)
	if limited < 0 {
		limited += 1.0
	}

	return limited
}

func (s *spa) limitMinutes(minutes float64) float64 {
	limited := minutes
	if limited < -20.0 {
		limited += 1440.0
	} else if limited > 20.0 {
		limited -= 1440.0
	}
	return limited
}

func (s *spa) dayfracToLocalHr(dayfrac float64, timezone float64) float64 {
	return 24.0 * s.limitZero2one(dayfrac+timezone/24.0)
}

func (s *spa) thirdOrderPolynomial(a float64, b float64, c float64, d float64, x float64) float64 {
	return ((a*x+b)*x+c)*x + d
}

func (s *spa) julianDay(year int, month int, day int, hour int, minute int, second float64, dut1 float64, tz float64) float64 {
	var dayDecimal, julianDay, a float64
	dayDecimal = float64(day) + (float64(hour)-tz+(float64(minute)+(second+dut1)/60.0)/60.0)/24.0

	if month < 3 {
		month += 12
		year--
	}

	julianDay = float64(int(365.25*(float64(year)+4716.0))) + float64(int(30.6001*(float64(month)+1))) + dayDecimal - 1524.5

	if julianDay > 2299160.0 {
		a = float64(year / 100)
		julianDay += 2 - a + float64(int(a/4))
	}

	return julianDay
}
func (s *spa) julianCentury(jd float64) float64 {
	return (jd - 2451545.0) / 36525.0
}

func (s *spa) julianEphemerisDay(jd float64, deltaT float64) float64 {
	return jd + deltaT/86400.0
}

func (s *spa) julianEphemerisCentury(jde float64) float64 {
	return (jde - 2451545.0) / 36525.0
}

func (s *spa) julianEphemerisMillennium(jce float64) float64 {
	return jce / 10.0
}

func (s *spa) earthPeriodicTermSummation(terms [][]float64, count int, jme float64) (sum float64) {
	sum = 0
	for i := 0; i < count; i++ {
		sum += terms[i][TermA] * math.Cos(terms[i][TermB]+terms[i][TermC]*jme)
	}
	return sum
}

func (s *spa) earthValues(termSum []float64, count int64, jme float64) (sum float64) {
	sum = 0
	for i := 0; i < int(count); i++ {
		sum += termSum[i] * math.Pow(jme, float64(i))
	}
	sum /= 1.0e8

	return sum
}

func (s *spa) earthHeliocentricLongitude(jme float64) float64 {
	sum := make([]float64, LCount)

	for i := 0; i < int(LCount); i++ {

		sum[i] = s.earthPeriodicTermSummation(LTerms[i], int(lSubcount[i]), jme)
	}
	return s.limitDegrees(s.rad2deg(s.earthValues(sum, LCount, jme)))

}

func (s *spa) earthHeliocentricLatitude(jme float64) float64 {

	sum := make([]float64, BCount)
	for i := 0; i < int(BCount); i++ {

		sum[i] = s.earthPeriodicTermSummation(BTerms[i], int(bSubcount[i]), jme)
	}

	return s.rad2deg(s.earthValues(sum, BCount, jme))

}

func (s *spa) earthRadiusVector(jme float64) float64 {
	sum := make([]float64, RCount)
	for i := 0; i < int(RCount); i++ {

		sum[i] = s.earthPeriodicTermSummation(RTerms[i], int(rSubcount[i]), jme)
	}
	return s.earthValues(sum, RCount, jme)

}

func (s *spa) geocentricLongitude(l float64) float64 {

	theta := l + 180.0

	if theta >= 360.0 {
		theta -= 360.0
	}

	return theta
}

func (s *spa) geocentricLatitude(b float64) float64 {
	return -b
}

func (s *spa) meanElongationMoonSun(jce float64) float64 {
	return s.thirdOrderPolynomial(1.0/189474.0, -0.0019142, 445267.11148, 297.85036, jce)
}

func (s *spa) meanAnomalySun(jce float64) float64 {
	return s.thirdOrderPolynomial(-1.0/300000.0, -0.0001603, 35999.05034, 357.52772, jce)
}

func (s *spa) meanAnomalyMoon(jce float64) float64 {
	return s.thirdOrderPolynomial(1.0/56250.0, 0.0086972, 477198.867398, 134.96298, jce)
}

func (s *spa) argumentLatitudeMoon(jce float64) float64 {
	return s.thirdOrderPolynomial(1.0/327270.0, -0.0036825, 483202.017538, 93.27191, jce)
}

func (s *spa) ascendingLongitudeMoon(jce float64) float64 {
	return s.thirdOrderPolynomial(1.0/450000.0, 0.0020708, -1934.136261, 125.04452, jce)
}

func (s *spa) xyTermSummation(i int, x []float64) (sum float64) {

	sum = 0
	for j := 0; j < int(TermYCount); j++ {

		sum += x[j] * float64(YTerms[i][j])
	}
	return sum
}

func (s *spa) nutationLongitudeAndObliquity(jce float64, x []float64) {
	var xyTermSum float64
	sumPsi := 0.
	sumEpsilon := 0.
	for i := 0; i < int(YCount); i++ {

		xyTermSum = s.deg2rad(s.xyTermSummation(i, x))
		sumPsi += (PeTerms[i][TermPsiA] + jce*PeTerms[i][TermPsiB]) * math.Sin(xyTermSum)
		sumEpsilon += (PeTerms[i][TermEpsC] + jce*PeTerms[i][TermEpsD]) * math.Cos(xyTermSum)
	}
	s.delPsi = sumPsi / 36000000.0
	s.delEpsilon = sumEpsilon / 36000000.0
}

func (s *spa) eclipticMeanObliquity(jme float64) float64 {
	u := jme / 10.0

	return 84381.448 + u*(-4680.93+u*(-1.55+u*(1999.25+u*(-51.38+u*(-249.67+
		u*(-39.05+u*(7.12+u*(27.87+u*(5.79+u*2.45)))))))))
}

func (s *spa) eclipticTrueObliquity(deltaEpsilon float64, epsilon0 float64) float64 {
	return deltaEpsilon + epsilon0/3600.0
}

func (s *spa) aberrationCorrection(r float64) float64 {
	return -20.4898 / (3600.0 * r)
}

func (s *spa) apparentSunLongitude(theta float64, deltaPsi float64, deltaTau float64) float64 {
	return theta + deltaPsi + deltaTau
}

func (s *spa) greenwichMeanSiderealTime(jd float64, jc float64) float64 {
	return s.limitDegrees(280.46061837 + 360.98564736629*(jd-2451545.0) +
		jc*jc*(0.000387933-jc/38710000.0))
}

func (s *spa) greenwichSiderealTime(nu0 float64, deltaPsi float64, epsilon float64) float64 {
	return nu0 + deltaPsi*math.Cos(s.deg2rad(epsilon))
}

func (s *spa) geocentricRightAscension(lamda float64, epsilon float64, beta float64) float64 {
	lamdaRad := s.deg2rad(lamda)
	epsilonRad := s.deg2rad(epsilon)

	return s.limitDegrees(s.rad2deg(math.Atan2(math.Sin(lamdaRad)*math.Cos(epsilonRad)-
		math.Tan(s.deg2rad(beta))*math.Sin(epsilonRad), math.Cos(lamdaRad))))
}

func (s *spa) geocentricDeclination(beta float64, epsilon float64, lambda float64) float64 {
	betaRad := s.deg2rad(beta)
	epsilonRad := s.deg2rad(epsilon)

	return s.rad2deg(math.Asin(math.Sin(betaRad)*math.Cos(epsilonRad) +
		math.Cos(betaRad)*math.Sin(epsilonRad)*math.Sin(s.deg2rad(lambda))))
}

func (s *spa) observerHourAngle(nu float64, longitude float64, alphaDeg float64) float64 {
	return s.limitDegrees(nu + longitude - alphaDeg)
}

func (s *spa) sunEquatorialHorizontalParallax(r float64) float64 {
	return 8.794 / (3600.0 * r)
}
func (s *spa) rightAscensionParallaxAndTopocentricDec(latitude float64, elevation float64, xi float64, h float64, delta float64) {
	var deltaAlphaRad float64
	latRad := s.deg2rad(latitude)
	xiRad := s.deg2rad(xi)
	hRad := s.deg2rad(h)
	deltaRad := s.deg2rad(delta)
	u := math.Atan(0.99664719 * math.Tan(latRad))
	y := 0.99664719*math.Sin(u) + elevation*math.Sin(latRad)/6378140.0
	x := math.Cos(u) + elevation*math.Cos(latRad)/6378140.0

	deltaAlphaRad = math.Atan2(-x*math.Sin(xiRad)*math.Sin(hRad), math.Cos(deltaRad)-x*math.Sin(xiRad)*math.Cos(hRad))

	s.deltaPrime = s.rad2deg(math.Atan2((math.Sin(deltaRad)-y*math.Sin(xiRad))*math.Cos(deltaAlphaRad),
		math.Cos(deltaRad)-x*math.Sin(xiRad)*math.Cos(hRad)))

	s.delAlpha = s.rad2deg(deltaAlphaRad)
}

func (s *spa) topocentricRightAscension(alphaDeg float64, deltaAlpha float64) float64 {
	return alphaDeg + deltaAlpha
}

func (s *spa) topocentricLocalHourAngle(h float64, deltaAlpha float64) float64 {
	return h - deltaAlpha
}

func (s *spa) topocentricElevationAngle(latitude float64, deltaPrime float64, hPrime float64) float64 {
	latRad := s.deg2rad(latitude)
	deltaPrimeRad := s.deg2rad(deltaPrime)

	return s.rad2deg(math.Asin(math.Sin(latRad)*math.Sin(deltaPrimeRad) +
		math.Cos(latRad)*math.Cos(deltaPrimeRad)*math.Cos(s.deg2rad(hPrime))))
}

func (s *spa) atmosphericRefractionCorrection(pressure float64, temperature float64, atmosRefract float64, e0 float64) float64 {
	delE := 0.

	if e0 >= -1*(SunRadius+atmosRefract) {
		delE = (pressure / 1010.0) * (283.0 / (273.0 + temperature)) * 1.02 / (60.0 * math.Tan(s.deg2rad(e0+10.3/(e0+5.11))))
	}
	return delE
}

func (s *spa) topocentricElevationAngleCorrected(e0 float64, deltaE float64) float64 {
	return e0 + deltaE
}

func (s *spa) topocentricZenithAngle(e float64) float64 {
	return 90.0 - e
}

func (s *spa) topocentricAzimuthAngleAstro(hPrime float64, latitude float64, deltaPrime float64) float64 {
	hPrimeRad := s.deg2rad(hPrime)
	latRad := s.deg2rad(latitude)

	return s.limitDegrees(s.rad2deg(math.Atan2(math.Sin(hPrimeRad),
		math.Cos(hPrimeRad)*math.Sin(latRad)-math.Tan(s.deg2rad(deltaPrime))*math.Cos(latRad))))
}

func (s *spa) topocentricAzimuthAngle(azimuthAstro float64) float64 {
	return s.limitDegrees(azimuthAstro + 180.0)
}

func (s *spa) surfaceIncidenceAngle(zenith float64, azimuthAstro float64, azmRotation float64, slope float64) float64 {
	zenithRad := s.deg2rad(zenith)
	slopeRad := s.deg2rad(slope)

	return s.rad2deg(math.Acos(math.Cos(zenithRad)*math.Cos(slopeRad) +
		math.Sin(slopeRad)*math.Sin(zenithRad)*math.Cos(s.deg2rad(azimuthAstro-azmRotation))))
}

func (s *spa) sunMeanLongitude(jme float64) float64 {
	return s.limitDegrees(280.4664567 + jme*(360007.6982779+jme*(0.03032028+
		jme*(1/49931.0+jme*(-1/15300.0+jme*(-1/2000000.0))))))
}

func (s *spa) eotf(m float64, alpha float64, delPsi float64, epsilon float64) float64 {
	return s.limitMinutes(4.0 * (m - 0.0057183 - alpha + delPsi*math.Cos(s.deg2rad(epsilon))))
}

func (s *spa) approxSunTransitTime(alphaZero float64, longitude float64, nu float64) float64 {
	return (alphaZero - longitude - nu) / 360.0
}

func (s *spa) sunHourAngleAtRiseSet(latitude float64, deltaZero float64, h0Prime float64) float64 {
	h0 := -99999.
	latitudeRad := s.deg2rad(latitude)
	deltaZeroRad := s.deg2rad(deltaZero)
	argument := (math.Sin(s.deg2rad(h0Prime)) - math.Sin(latitudeRad)*math.Sin(deltaZeroRad)) /
		(math.Cos(latitudeRad) * math.Cos(deltaZeroRad))

	if math.Abs(argument) <= 1 {
		h0 = s.limitDegrees180(s.rad2deg(math.Acos(argument)))
	}

	return h0
}

func (s *spa) approxSunRiseAndSet(h0 float64) {
	h0Dfrac := h0 / 360.0
	s.mRts[SunRise] = s.limitZero2one(s.mRts[SunTransit] - h0Dfrac)
	s.mRts[SunSet] = s.limitZero2one(s.mRts[SunTransit] + h0Dfrac)
	s.mRts[SunTransit] = s.limitZero2one(s.mRts[SunTransit])
}

func (s *spa) rtsAlphaDeltaPrime(ad []float64, n float64) float64 {
	a := ad[JdZero] - ad[JdMinus]
	b := ad[JdPlus] - ad[JdZero]

	if math.Abs(a) >= 2.0 {
		a = s.limitZero2one(a)
	}
	if math.Abs(b) >= 2.0 {
		b = s.limitZero2one(b)
	}

	return ad[JdZero] + n*(a+b+(b-a)*n)/2.0
}

func (s *spa) rtsSunAltitude(latitude float64, deltaPrime float64, hPrime float64) float64 {
	latitudeRad := s.deg2rad(latitude)
	deltaPrimeRad := s.deg2rad(deltaPrime)

	return s.rad2deg(math.Asin(math.Sin(latitudeRad)*math.Sin(deltaPrimeRad) +
		math.Cos(latitudeRad)*math.Cos(deltaPrimeRad)*math.Cos(s.deg2rad(hPrime))))
}

func (s *spa) sunRiseAndSet(mRts []float64, hRts []float64, deltaPrime []float64, latitude float64, hPrime []float64, h0Prime float64, sun int) float64 {
	return mRts[sun] + (hRts[sun]-h0Prime)/
		(360.0*math.Cos(s.deg2rad(deltaPrime[sun]))*math.Cos(s.deg2rad(latitude))*math.Sin(s.deg2rad(hPrime[sun])))
}

////////////////////////////////////////////////////////////////////////////////////////////////
// Calculate required SPA parameters to get the right ascension (alpha) and declination (delta)
// Note: JD must be already calculated and in structure
////////////////////////////////////////////////////////////////////////////////////////////////
func (s *spa) calculateGeocentricSunRightAscensionAndDeclination() {
	x := make([]float64, TermXCount)

	s.jc = s.julianCentury(s.jd)

	s.jde = s.julianEphemerisDay(s.jd, s.deltaT)
	s.jce = s.julianEphemerisCentury(s.jde)
	s.jme = s.julianEphemerisMillennium(s.jce)

	s.l = s.earthHeliocentricLongitude(s.jme)
	s.b = s.earthHeliocentricLatitude(s.jme)
	s.r = s.earthRadiusVector(s.jme)

	s.theta = s.geocentricLongitude(s.l)
	s.beta = s.geocentricLatitude(s.b)

	x[TermX0], s.x0 = s.meanElongationMoonSun(s.jce), s.meanElongationMoonSun(s.jce)
	x[TermX1], s.x1 = s.meanAnomalySun(s.jce), s.meanAnomalySun(s.jce)
	x[TermX2], s.x2 = s.meanAnomalyMoon(s.jce), s.meanAnomalyMoon(s.jce)
	x[TermX3], s.x3 = s.argumentLatitudeMoon(s.jce), s.argumentLatitudeMoon(s.jce)
	x[TermX4], s.x4 = s.ascendingLongitudeMoon(s.jce), s.ascendingLongitudeMoon(s.jce)

	s.nutationLongitudeAndObliquity(s.jce, x)

	s.epsilon0 = s.eclipticMeanObliquity(s.jme)
	s.epsilon = s.eclipticTrueObliquity(s.delEpsilon, s.epsilon0)

	s.delTau = s.aberrationCorrection(s.r)
	s.lamda = s.apparentSunLongitude(s.theta, s.delPsi, s.delTau)
	s.nu0 = s.greenwichMeanSiderealTime(s.jd, s.jc)
	s.nu = s.greenwichSiderealTime(s.nu0, s.delPsi, s.epsilon)

	s.alpha = s.geocentricRightAscension(s.lamda, s.epsilon, s.beta)
	s.delta = s.geocentricDeclination(s.beta, s.epsilon, s.lamda)
}

////////////////////////////////////////////////////////////////////////
// Calculate Equation of Time (EOT) and Sun Rise, Transit, & Set (RTS)
////////////////////////////////////////////////////////////////////////

func (s *spa) calculateEotAndSunRiseTransitSet() {
	//spa_data
	// sun_rts
	var nu, m, h0, n float64
	alpha := make([]float64, JdCount)
	delta := make([]float64, JdCount)

	nuRts := make([]float64, SunCount)
	hRts := make([]float64, SunCount)
	alphaPrime := make([]float64, SunCount)
	deltaPrime := make([]float64, SunCount)
	hPrime := make([]float64, SunCount)
	h0Prime := -1 * (SunRadius + s.atmosRefract)

	sunRts := *s
	sunRts.mRts = make([]float64, SunCount)
	m = s.sunMeanLongitude(s.jme)
	s.eot = s.eotf(m, s.alpha, s.delPsi, s.epsilon)

	sunRts.hour = 0
	sunRts.minute = 0
	sunRts.second = 0
	sunRts.deltaUt1 = 0
	sunRts.timezone = 0.0

	sunRts.jd = s.julianDay(sunRts.year, sunRts.month, sunRts.day, sunRts.hour,
		sunRts.minute, sunRts.second, sunRts.deltaUt1, sunRts.timezone)

	sunRts.calculateGeocentricSunRightAscensionAndDeclination()
	nu = sunRts.nu

	sunRts.deltaT = 0
	sunRts.jd--
	for i := 0; i < JdCount; i++ {
		sunRts.calculateGeocentricSunRightAscensionAndDeclination()
		alpha[i] = sunRts.alpha
		delta[i] = sunRts.delta
		sunRts.jd++
	}

	sunRts.mRts[SunTransit] = s.approxSunTransitTime(alpha[JdZero], s.longitude, nu)
	h0 = s.sunHourAngleAtRiseSet(s.latitude, delta[JdZero], h0Prime)

	if h0 >= 0 {

		sunRts.approxSunRiseAndSet(h0)
		for i := 0; i < SunCount; i++ {

			nuRts[i] = nu + 360.985647*sunRts.mRts[i]

			n = sunRts.mRts[i] + sunRts.deltaT/86400.0
			alphaPrime[i] = sunRts.rtsAlphaDeltaPrime(alpha, n)
			deltaPrime[i] = sunRts.rtsAlphaDeltaPrime(delta, n)

			hPrime[i] = s.limitDegrees180pm(nuRts[i] + s.longitude - alphaPrime[i])

			hRts[i] = s.rtsSunAltitude(s.latitude, deltaPrime[i], hPrime[i])
		}

		s.srha = hPrime[SunRise]
		s.ssha = hPrime[SunSet]
		s.sta = hRts[SunTransit]

		s.suntransit = s.dayfracToLocalHr(sunRts.mRts[SunTransit]-hPrime[SunTransit]/360.0,
			s.timezone)

		s.sunrise = s.dayfracToLocalHr(s.sunRiseAndSet(sunRts.mRts, hRts, deltaPrime,
			s.latitude, hPrime, h0Prime, SunRise), s.timezone)

		s.sunset = s.dayfracToLocalHr(s.sunRiseAndSet(sunRts.mRts, hRts, deltaPrime,
			s.latitude, hPrime, h0Prime, SunSet), s.timezone)

	} else {
		s.srha, s.ssha, s.sta, s.suntransit, s.sunrise, s.sunset = -99999, -99999, -99999, -99999, -99999, -99999
	}

}

func (s *spa) validate() error {
	if (s.year < -2000) || (s.year > 6000) {
		return errors.New("invalid year")
	}
	if (s.month < 1) || (s.month > 12) {
		return errors.New("invalid month")
	}
	if (s.day < 1) || (s.day > 31) {
		return errors.New("invalid day")
	}
	if (s.hour < 0) || (s.hour > 24) {
		return errors.New("invalid hour")
	}
	if (s.minute < 0) || (s.minute > 59) {
		return errors.New("invalid minute")
	}
	if (s.second < 0) || (s.second >= 60) {
		return errors.New("invalid second")
	}
	if (s.pressure < 0) || (s.pressure > 5000) {
		return errors.New("invalid pressure")
	}
	if (s.temperature <= -273) || (s.temperature > 6000) {
		return errors.New("invalid temperature")
	}
	if (s.deltaUt1 <= -1) || (s.deltaUt1 >= 1) {
		return errors.New("invalid UTC / UT difference (deltaUt1)")
	}
	if (s.hour == 24) && (s.minute > 0) {
		return errors.New("invalid hour/minute")
	}
	if (s.hour == 24) && (s.second > 0) {
		return errors.New("invalid hour/second")
	}

	if math.Abs(s.deltaT) > 8000 {
		return errors.New("invalid difference between earth rotation time and terrestrial time (deltaT)")
	}
	if math.Abs(s.timezone) > 18 {
		return errors.New("invalid timezone")
	}
	if math.Abs(s.longitude) > 180 {
		return errors.New("invalid longitude")
	}
	if math.Abs(s.latitude) > 90 {
		return errors.New("invalid latitude")
	}
	if math.Abs(s.atmosRefract) > 5 {
		return errors.New("invalid atmospheric refraction (atmosRefract)")
	}
	if s.elevation < -6500000 {
		return errors.New("invalid elevation")
	}

	if (s.function == SpaZaInc) || (s.function == SpaAll) {
		if math.Abs(s.slope) > 360 {
			return errors.New("invalid surface slope")
		}
		if math.Abs(s.azmRotation) > 360 {
			return errors.New("invalid surface azimuth rotation")
		}
	}

	return nil
}
