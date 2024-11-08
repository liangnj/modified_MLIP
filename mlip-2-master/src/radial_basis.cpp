/*   This software is called MLIP for Machine Learning Interatomic Potentials.
 *   MLIP can only be used for non-commercial research and cannot be re-distributed.
 *   The use of MLIP must be acknowledged by citing approriate references.
 *   See the LICENSE file for details.
 *
 *   This file contributors: Alexander Shapeev
 */

#include "radial_basis.h"


#include <cmath>

using namespace std;


void AnyRadialBasis::ReadRadialBasis(std::ifstream & ifs)
{
	if ((!ifs.is_open()) || (ifs.eof()))
		ERROR("RadialBasis::ReadRadialBasis: Can't load radial basis");

	string tmpstr;

	// reading min_dist / scaling
	ifs >> tmpstr;
	if (tmpstr == "scaling") {
		ifs.ignore(2);
		ifs >> scaling;
		if (ifs.fail())
			ERROR("Error reading .mtp file");
		ifs >> tmpstr;
	}

	if (tmpstr != "min_dist")
		ERROR("Error reading .mtp file");
	ifs.ignore(2);
	ifs >> min_dist;
	if (ifs.fail())
		ERROR("Error reading .mtp file");

	// reading max_dist 
	ifs >> tmpstr;
	if (tmpstr != "max_dist")
		ERROR("Error reading .mtp file");
	ifs.ignore(2);
	ifs >> max_dist;
	if (ifs.fail())
		ERROR("Error reading .mtp file");

	// reading rb_size 
	ifs >> tmpstr;
	if (tmpstr != "radial_basis_size")
		ERROR("Error reading .mtp file");
	ifs.ignore(2);
	ifs >> rb_size;
	if (ifs.fail())
		ERROR("Error reading .mtp file");

	rb_vals.resize(rb_size);
	rb_ders.resize(rb_size);
}

void AnyRadialBasis::WriteRadialBasis(std::ofstream & ofs)
{
	if (!ofs.is_open())
		ERROR("RadialBasis::WriteRadialBasis: Output stream isn't open");

	ofs << "radial_basis_type = " << GetRBTypeString() << '\n';
	if (scaling != 1.0)
		ofs << "\tscaling = " << scaling << '\n';
	ofs << "\tmin_dist = " << min_dist << '\n'
		<< "\tmax_dist = " << max_dist << '\n'
		<< "\tradial_basis_size = " << rb_size << '\n';
}

AnyRadialBasis::AnyRadialBasis(double _min_dist, double _max_dist, int _size)
	: rb_size(_size), min_dist(_min_dist), max_dist(_max_dist)
{
	rb_vals.resize(rb_size);
	rb_ders.resize(rb_size);
}

AnyRadialBasis::AnyRadialBasis(std::ifstream & ifs)
{
	ReadRadialBasis(ifs);
}

void RadialBasis_Shapeev::InitShapeevRB()
{
#ifdef MLIP_DEBUG
	if (rb_size > ALLOCATED_DEGREE) {
		ERROR("error: RadialBasis::MAX_DEGREE > RadialBasis::ALLOCATED_DEGREE");
	}
#endif

	rb_coeffs[0][0] = (sqrt(21)*pow(max_dist, 2)*pow(min_dist, 2)) / (2.*pow(max_dist - min_dist, 2));
	rb_coeffs[0][1] = -((sqrt(21)*max_dist*pow(min_dist, 2)) / pow(max_dist - min_dist, 2));
	rb_coeffs[0][2] = (sqrt(21)*pow(min_dist, 2)) / (2.*pow(max_dist - min_dist, 2));
	rb_coeffs[0][3] = 0;
	rb_coeffs[0][4] = 0;
	rb_coeffs[0][5] = 0;
	rb_coeffs[0][6] = 0;
	rb_coeffs[0][7] = 0;
	rb_coeffs[0][8] = 0;
	rb_coeffs[0][9] = 0;
	rb_coeffs[0][10] = 0;
	rb_coeffs[0][11] = 0;
	rb_coeffs[0][12] = 0;
	rb_coeffs[0][13] = 0;
	rb_coeffs[1][0] = (-3 * sqrt(7)*pow(max_dist, 2)*pow(min_dist, 2)*(max_dist + 3 * min_dist)) / (4.*pow(max_dist - min_dist, 3));
	rb_coeffs[1][1] = (9 * sqrt(7)*max_dist*pow(min_dist, 2)*(max_dist + min_dist)) / (2.*pow(max_dist - min_dist, 3));
	rb_coeffs[1][2] = (9 * sqrt(7)*pow(min_dist, 2)*(3 * max_dist + min_dist)) / (4.*pow(-max_dist + min_dist, 3));
	rb_coeffs[1][3] = (3 * sqrt(7)*pow(min_dist, 2)) / pow(max_dist - min_dist, 3);
	rb_coeffs[1][4] = 0;
	rb_coeffs[1][5] = 0;
	rb_coeffs[1][6] = 0;
	rb_coeffs[1][7] = 0;
	rb_coeffs[1][8] = 0;
	rb_coeffs[1][9] = 0;
	rb_coeffs[1][10] = 0;
	rb_coeffs[1][11] = 0;
	rb_coeffs[1][12] = 0;
	rb_coeffs[1][13] = 0;
	rb_coeffs[2][0] = (sqrt(3.6666666666666665)*pow(max_dist, 2)*pow(min_dist, 2)*(pow(max_dist, 2) + 7 * max_dist*min_dist + 7 * pow(min_dist, 2))) / pow(max_dist - min_dist, 4);
	rb_coeffs[2][1] = -((sqrt(3.6666666666666665)*max_dist*pow(min_dist, 2)*(11 * pow(max_dist, 2) + 35 * max_dist*min_dist + 14 * pow(min_dist, 2))) / pow(max_dist - min_dist, 4));
	rb_coeffs[2][2] = (sqrt(3.6666666666666665)*pow(min_dist, 2)*(34 * pow(max_dist, 2) + 49 * max_dist*min_dist + 7 * pow(min_dist, 2))) / pow(max_dist - min_dist, 4);
	rb_coeffs[2][3] = -((sqrt(33)*pow(min_dist, 2)*(13 * max_dist + 7 * min_dist)) / pow(max_dist - min_dist, 4));
	rb_coeffs[2][4] = (5 * sqrt(33)*pow(min_dist, 2)) / pow(max_dist - min_dist, 4);
	rb_coeffs[2][5] = 0;
	rb_coeffs[2][6] = 0;
	rb_coeffs[2][7] = 0;
	rb_coeffs[2][8] = 0;
	rb_coeffs[2][9] = 0;
	rb_coeffs[2][10] = 0;
	rb_coeffs[2][11] = 0;
	rb_coeffs[2][12] = 0;
	rb_coeffs[2][13] = 0;
	rb_coeffs[3][0] = (-3 * sqrt(6.5)*pow(max_dist, 2)*pow(min_dist, 2)*(pow(max_dist, 3) + 12 * pow(max_dist, 2)*min_dist + 28 * max_dist*pow(min_dist, 2) + 14 * pow(min_dist, 3))) / (4.*pow(max_dist - min_dist, 5));
	rb_coeffs[3][1] = (3 * sqrt(6.5)*max_dist*pow(min_dist, 2)*(17 * pow(max_dist, 3) + 104 * pow(max_dist, 2)*min_dist + 126 * max_dist*pow(min_dist, 2) + 28 * pow(min_dist, 3))) / (4.*pow(max_dist - min_dist, 5));
	rb_coeffs[3][2] = (3 * sqrt(6.5)*pow(min_dist, 2)*(43 * pow(max_dist, 3) + 141 * pow(max_dist, 2)*min_dist + 84 * max_dist*pow(min_dist, 2) + 7 * pow(min_dist, 3))) / (2.*pow(-max_dist + min_dist, 5));
	rb_coeffs[3][3] = (15 * sqrt(6.5)*pow(min_dist, 2)*(18 * pow(max_dist, 2) + 30 * max_dist*min_dist + 7 * pow(min_dist, 2))) / (2.*pow(max_dist - min_dist, 5));
	rb_coeffs[3][4] = (165 * sqrt(6.5)*pow(min_dist, 2)*(3 * max_dist + 2 * min_dist)) / (4.*pow(-max_dist + min_dist, 5));
	rb_coeffs[3][5] = (165 * sqrt(6.5)*pow(min_dist, 2)) / (4.*pow(max_dist - min_dist, 5));
	rb_coeffs[3][6] = 0;
	rb_coeffs[3][7] = 0;
	rb_coeffs[3][8] = 0;
	rb_coeffs[3][9] = 0;
	rb_coeffs[3][10] = 0;
	rb_coeffs[3][11] = 0;
	rb_coeffs[3][12] = 0;
	rb_coeffs[3][13] = 0;
	rb_coeffs[4][0] = (sqrt(0.6)*pow(max_dist, 2)*pow(min_dist, 2)*(5 * pow(max_dist, 4) + 90 * pow(max_dist, 3)*min_dist + 360 * pow(max_dist, 2)*pow(min_dist, 2) + 420 * max_dist*pow(min_dist, 3) + 126 * pow(min_dist, 4))) / (2.*pow(max_dist - min_dist, 6));
	rb_coeffs[4][1] = (-3 * sqrt(0.6)*max_dist*pow(min_dist, 2)*(20 * pow(max_dist, 4) + 195 * pow(max_dist, 3)*min_dist + 450 * pow(max_dist, 2)*pow(min_dist, 2) + 294 * max_dist*pow(min_dist, 3) + 42 * pow(min_dist, 4))) / pow(max_dist - min_dist, 6);
	rb_coeffs[4][2] = (3 * sqrt(0.6)*pow(min_dist, 2)*(295 * pow(max_dist, 4) + 1680 * pow(max_dist, 3)*min_dist + 2232 * pow(max_dist, 2)*pow(min_dist, 2) + 756 * max_dist*pow(min_dist, 3) + 42 * pow(min_dist, 4))) / (2.*pow(max_dist - min_dist, 6));
	rb_coeffs[4][3] = (-22 * sqrt(0.6)*pow(min_dist, 2)*(65 * pow(max_dist, 3) + 216 * pow(max_dist, 2)*min_dist + 153 * max_dist*pow(min_dist, 2) + 21 * pow(min_dist, 3))) / pow(max_dist - min_dist, 6);
	rb_coeffs[4][4] = (33 * sqrt(0.6)*pow(min_dist, 2)*(137 * pow(max_dist, 2) + 246 * max_dist*min_dist + 72 * pow(min_dist, 2))) / (2.*pow(max_dist - min_dist, 6));
	rb_coeffs[4][5] = (-429 * sqrt(0.6)*pow(min_dist, 2)*(4 * max_dist + 3 * min_dist)) / pow(max_dist - min_dist, 6);
	rb_coeffs[4][6] = (1001 * sqrt(0.6)*pow(min_dist, 2)) / (2.*pow(max_dist - min_dist, 6));
	rb_coeffs[4][7] = 0;
	rb_coeffs[4][8] = 0;
	rb_coeffs[4][9] = 0;
	rb_coeffs[4][10] = 0;
	rb_coeffs[4][11] = 0;
	rb_coeffs[4][12] = 0;
	rb_coeffs[4][13] = 0;
	rb_coeffs[5][0] = -(sqrt(62.333333333333336)*pow(max_dist, 2)*pow(min_dist, 2)*(pow(max_dist, 5) + 25 * pow(max_dist, 4)*min_dist + 150 * pow(max_dist, 3)*pow(min_dist, 2) + 300 * pow(max_dist, 2)*pow(min_dist, 3) + 210 * max_dist*pow(min_dist, 4) + 42 * pow(min_dist, 5))) / (4.*pow(max_dist - min_dist, 7));
	rb_coeffs[5][1] = (sqrt(62.333333333333336)*max_dist*pow(min_dist, 2)*(16 * pow(max_dist, 5) + 225 * pow(max_dist, 4)*min_dist + 825 * pow(max_dist, 3)*pow(min_dist, 2) + 1020 * pow(max_dist, 2)*pow(min_dist, 3) + 420 * max_dist*pow(min_dist, 4) + 42 * pow(min_dist, 5))) / (2.*pow(max_dist - min_dist, 7));
	rb_coeffs[5][2] = -(sqrt(561)*pow(min_dist, 2)*(107 * pow(max_dist, 5) + 925 * pow(max_dist, 4)*min_dist + 2120 * pow(max_dist, 3)*pow(min_dist, 2) + 1580 * pow(max_dist, 2)*pow(min_dist, 3) + 350 * max_dist*pow(min_dist, 4) + 14 * pow(min_dist, 5))) / (4.*pow(max_dist - min_dist, 7));
	rb_coeffs[5][3] = (5 * sqrt(62.333333333333336)*pow(min_dist, 2)*(73 * pow(max_dist, 4) + 397 * pow(max_dist, 3)*min_dist + 555 * pow(max_dist, 2)*pow(min_dist, 2) + 228 * max_dist*pow(min_dist, 3) + 21 * pow(min_dist, 4))) / pow(max_dist - min_dist, 7);
	rb_coeffs[5][4] = (-65 * sqrt(62.333333333333336)*pow(min_dist, 2)*(53 * pow(max_dist, 3) + 177 * pow(max_dist, 2)*min_dist + 138 * max_dist*pow(min_dist, 2) + 24 * pow(min_dist, 3))) / (4.*pow(max_dist - min_dist, 7));
	rb_coeffs[5][5] = (91 * sqrt(561)*pow(min_dist, 2)*(8 * pow(max_dist, 2) + 15 * max_dist*min_dist + 5 * pow(min_dist, 2))) / (2.*pow(max_dist - min_dist, 7));
	rb_coeffs[5][6] = (-91 * sqrt(62.333333333333336)*pow(min_dist, 2)*(31 * max_dist + 25 * min_dist)) / (4.*pow(max_dist - min_dist, 7));
	rb_coeffs[5][7] = (182 * sqrt(62.333333333333336)*pow(min_dist, 2)) / pow(max_dist - min_dist, 7);
	rb_coeffs[5][8] = 0;
	rb_coeffs[5][9] = 0;
	rb_coeffs[5][10] = 0;
	rb_coeffs[5][11] = 0;
	rb_coeffs[5][12] = 0;
	rb_coeffs[5][13] = 0;
	rb_coeffs[6][0] = (sqrt(4.071428571428571)*pow(max_dist, 2)*pow(min_dist, 2)*(pow(max_dist, 6) + 33 * pow(max_dist, 5)*min_dist + 275 * pow(max_dist, 4)*pow(min_dist, 2) + 825 * pow(max_dist, 3)*pow(min_dist, 3) + 990 * pow(max_dist, 2)*pow(min_dist, 4) + 462 * max_dist*pow(min_dist, 5) + 66 * pow(min_dist, 6))) / pow(max_dist - min_dist, 8);
	rb_coeffs[6][1] = -((sqrt(4.071428571428571)*max_dist*pow(min_dist, 2)*(41 * pow(max_dist, 6) + 781 * pow(max_dist, 5)*min_dist + 4125 * pow(max_dist, 4)*pow(min_dist, 2) + 8085 * pow(max_dist, 3)*pow(min_dist, 3) + 6270 * pow(max_dist, 2)*pow(min_dist, 4) + 1782 * max_dist*pow(min_dist, 5) + 132 * pow(min_dist, 6))) / pow(max_dist - min_dist, 8));
	rb_coeffs[6][2] = (3 * sqrt(16.285714285714285)*pow(min_dist, 2)*(89 * pow(max_dist, 6) + 1078 * pow(max_dist, 5)*min_dist + 3740 * pow(max_dist, 4)*pow(min_dist, 2) + 4785 * pow(max_dist, 3)*pow(min_dist, 3) + 2310 * pow(max_dist, 2)*pow(min_dist, 4) + 363 * max_dist*pow(min_dist, 5) + 11 * pow(min_dist, 6))) / pow(max_dist - min_dist, 8);
	rb_coeffs[6][3] = (-13 * sqrt(16.285714285714285)*pow(min_dist, 2)*(124 * pow(max_dist, 5) + 990 * pow(max_dist, 4)*min_dist + 2255 * pow(max_dist, 3)*pow(min_dist, 2) + 1815 * pow(max_dist, 2)*pow(min_dist, 3) + 495 * max_dist*pow(min_dist, 4) + 33 * pow(min_dist, 5))) / pow(max_dist - min_dist, 8);
	rb_coeffs[6][4] = (65 * sqrt(4.071428571428571)*pow(min_dist, 2)*(161 * pow(max_dist, 4) + 847 * pow(max_dist, 3)*min_dist + 1221 * pow(max_dist, 2)*pow(min_dist, 2) + 561 * max_dist*pow(min_dist, 3) + 66 * pow(min_dist, 4))) / pow(max_dist - min_dist, 8);
	rb_coeffs[6][5] = (-39 * sqrt(4.071428571428571)*pow(min_dist, 2)*(497 * pow(max_dist, 3) + 1661 * pow(max_dist, 2)*min_dist + 1375 * max_dist*pow(min_dist, 2) + 275 * pow(min_dist, 3))) / pow(max_dist - min_dist, 8);
	rb_coeffs[6][6] = (26 * sqrt(16.285714285714285)*pow(min_dist, 2)*(394 * pow(max_dist, 2) + 759 * max_dist*min_dist + 275 * pow(min_dist, 2))) / pow(max_dist - min_dist, 8);
	rb_coeffs[6][7] = (-442 * sqrt(16.285714285714285)*pow(min_dist, 2)*(13 * max_dist + 11 * min_dist)) / pow(max_dist - min_dist, 8);
	rb_coeffs[6][8] = (1326 * sqrt(16.285714285714285)*pow(min_dist, 2)) / pow(max_dist - min_dist, 8);
	rb_coeffs[6][9] = 0;
	rb_coeffs[6][10] = 0;
	rb_coeffs[6][11] = 0;
	rb_coeffs[6][12] = 0;
	rb_coeffs[6][13] = 0;
	rb_coeffs[7][0] = -(sqrt(273)*pow(max_dist, 2)*pow(min_dist, 2)*(pow(max_dist, 7) + 42 * pow(max_dist, 6)*min_dist + 462 * pow(max_dist, 5)*pow(min_dist, 2) + 1925 * pow(max_dist, 4)*pow(min_dist, 3) + 3465 * pow(max_dist, 3)*pow(min_dist, 4) + 2772 * pow(max_dist, 2)*pow(min_dist, 5) + 924 * max_dist*pow(min_dist, 6) + 99 * pow(min_dist, 7))) / (8.*pow(max_dist - min_dist, 9));
	rb_coeffs[7][1] = (3 * sqrt(273)*max_dist*pow(min_dist, 2)*(17 * pow(max_dist, 7) + 420 * pow(max_dist, 6)*min_dist + 3003 * pow(max_dist, 5)*pow(min_dist, 2) + 8470 * pow(max_dist, 4)*pow(min_dist, 3) + 10395 * pow(max_dist, 3)*pow(min_dist, 4) + 5544 * pow(max_dist, 2)*pow(min_dist, 5) + 1155 * max_dist*pow(min_dist, 6) + 66 * pow(min_dist, 7))) / (8.*pow(max_dist - min_dist, 9));
	rb_coeffs[7][2] = (-3 * sqrt(273)*pow(min_dist, 2)*(278 * pow(max_dist, 7) + 4473 * pow(max_dist, 6)*min_dist + 21714 * pow(max_dist, 5)*pow(min_dist, 2) + 41965 * pow(max_dist, 4)*pow(min_dist, 3) + 34650 * pow(max_dist, 3)*pow(min_dist, 4) + 11781 * pow(max_dist, 2)*pow(min_dist, 5) + 1386 * max_dist*pow(min_dist, 6) + 33 * pow(min_dist, 7))) / (8.*pow(max_dist - min_dist, 9));
	rb_coeffs[7][3] = (7 * sqrt(273)*pow(min_dist, 2)*(917 * pow(max_dist, 6) + 10038 * pow(max_dist, 5)*min_dist + 33495 * pow(max_dist, 4)*pow(min_dist, 2) + 43780 * pow(max_dist, 3)*pow(min_dist, 3) + 23265 * pow(max_dist, 2)*pow(min_dist, 4) + 4554 * max_dist*pow(min_dist, 5) + 231 * pow(min_dist, 6))) / (8.*pow(max_dist - min_dist, 9));
	rb_coeffs[7][4] = (-105 * sqrt(273)*pow(min_dist, 2)*(259 * pow(max_dist, 5) + 1953 * pow(max_dist, 4)*min_dist + 4422 * pow(max_dist, 3)*pow(min_dist, 2) + 3740 * pow(max_dist, 2)*pow(min_dist, 3) + 1155 * max_dist*pow(min_dist, 4) + 99 * pow(min_dist, 5))) / (8.*pow(max_dist - min_dist, 9));
	rb_coeffs[7][5] = (21 * sqrt(273)*pow(min_dist, 2)*(1624 * pow(max_dist, 4) + 8328 * pow(max_dist, 3)*min_dist + 12243 * pow(max_dist, 2)*pow(min_dist, 2) + 6050 * max_dist*pow(min_dist, 3) + 825 * pow(min_dist, 4))) / (4.*pow(max_dist - min_dist, 9));
	rb_coeffs[7][6] = (-119 * sqrt(273)*pow(min_dist, 2)*(436 * pow(max_dist, 3) + 1455 * pow(max_dist, 2)*min_dist + 1254 * max_dist*pow(min_dist, 2) + 275 * pow(min_dist, 3))) / (4.*pow(max_dist - min_dist, 9));
	rb_coeffs[7][7] = (153 * sqrt(273)*pow(min_dist, 2)*(307 * pow(max_dist, 2) + 602 * max_dist*min_dist + 231 * pow(min_dist, 2))) / (4.*pow(max_dist - min_dist, 9));
	rb_coeffs[7][8] = (-2907 * sqrt(273)*pow(min_dist, 2)*(8 * max_dist + 7 * min_dist)) / (4.*pow(max_dist - min_dist, 9));
	rb_coeffs[7][9] = (4845 * sqrt(273)*pow(min_dist, 2)) / (4.*pow(max_dist - min_dist, 9));
	rb_coeffs[7][10] = 0;
	rb_coeffs[7][11] = 0;
	rb_coeffs[7][12] = 0;
	rb_coeffs[7][13] = 0;
	rb_coeffs[8][0] = (sqrt(161)*pow(max_dist, 2)*pow(min_dist, 2)*(pow(max_dist, 8) + 52 * pow(max_dist, 7)*min_dist + 728 * pow(max_dist, 6)*pow(min_dist, 2) + 4004 * pow(max_dist, 5)*pow(min_dist, 3) + 10010 * pow(max_dist, 4)*pow(min_dist, 4) + 12012 * pow(max_dist, 3)*pow(min_dist, 5) + 6864 * pow(max_dist, 2)*pow(min_dist, 6) + 1716 * max_dist*pow(min_dist, 7) + 143 * pow(min_dist, 8))) / (6.*pow(max_dist - min_dist, 10));
	rb_coeffs[8][1] = -(sqrt(161)*max_dist*pow(min_dist, 2)*(31 * pow(max_dist, 8) + 962 * pow(max_dist, 7)*min_dist + 8918 * pow(max_dist, 6)*pow(min_dist, 2) + 34034 * pow(max_dist, 5)*pow(min_dist, 3) + 60060 * pow(max_dist, 4)*pow(min_dist, 4) + 50622 * pow(max_dist, 3)*pow(min_dist, 5) + 19734 * pow(max_dist, 2)*pow(min_dist, 6) + 3146 * max_dist*pow(min_dist, 7) + 143 * pow(min_dist, 8))) / (3.*pow(max_dist - min_dist, 10));
	rb_coeffs[8][2] = (sqrt(161)*pow(min_dist, 2)*(1241 * pow(max_dist, 8) + 25532 * pow(max_dist, 7)*min_dist + 164528 * pow(max_dist, 6)*pow(min_dist, 2) + 444444 * pow(max_dist, 5)*pow(min_dist, 3) + 553410 * pow(max_dist, 4)*pow(min_dist, 4) + 320892 * pow(max_dist, 3)*pow(min_dist, 5) + 81224 * pow(max_dist, 2)*pow(min_dist, 6) + 7436 * max_dist*pow(min_dist, 7) + 143 * pow(min_dist, 8))) / (6.*pow(max_dist - min_dist, 10));
	rb_coeffs[8][3] = (-10 * sqrt(161)*pow(min_dist, 2)*(591 * pow(max_dist, 7) + 8463 * pow(max_dist, 6)*min_dist + 38675 * pow(max_dist, 5)*pow(min_dist, 2) + 73931 * pow(max_dist, 4)*pow(min_dist, 3) + 63635 * pow(max_dist, 3)*pow(min_dist, 4) + 24167 * pow(max_dist, 2)*pow(min_dist, 5) + 3575 * max_dist*pow(min_dist, 6) + 143 * pow(min_dist, 7))) / (3.*pow(max_dist - min_dist, 10));
	rb_coeffs[8][4] = (10 * sqrt(161)*pow(min_dist, 2)*(3150 * pow(max_dist, 6) + 32032 * pow(max_dist, 5)*min_dist + 103792 * pow(max_dist, 4)*pow(min_dist, 2) + 137566 * pow(max_dist, 3)*pow(min_dist, 3) + 77935 * pow(max_dist, 2)*pow(min_dist, 4) + 17446 * max_dist*pow(min_dist, 5) + 1144 * pow(min_dist, 6))) / (3.*pow(max_dist - min_dist, 10));
	rb_coeffs[8][5] = (-68 * sqrt(161)*pow(min_dist, 2)*(1498 * pow(max_dist, 5) + 10816 * pow(max_dist, 4)*min_dist + 24349 * pow(max_dist, 3)*pow(min_dist, 2) + 21307 * pow(max_dist, 2)*pow(min_dist, 3) + 7150 * max_dist*pow(min_dist, 4) + 715 * pow(min_dist, 5))) / (3.*pow(max_dist - min_dist, 10));
	rb_coeffs[8][6] = (34 * sqrt(161)*pow(min_dist, 2)*(6102 * pow(max_dist, 4) + 30654 * pow(max_dist, 3)*min_dist + 45656 * pow(max_dist, 2)*pow(min_dist, 2) + 23738 * max_dist*pow(min_dist, 3) + 3575 * pow(min_dist, 4))) / (3.*pow(max_dist - min_dist, 10));
	rb_coeffs[8][7] = (-1292 * sqrt(161)*pow(min_dist, 2)*(207 * pow(max_dist, 3) + 689 * pow(max_dist, 2)*min_dist + 611 * max_dist*pow(min_dist, 2) + 143 * pow(min_dist, 3))) / (3.*pow(max_dist - min_dist, 10));
	rb_coeffs[8][8] = (1615 * sqrt(161)*pow(min_dist, 2)*(131 * pow(max_dist, 2) + 260 * max_dist*min_dist + 104 * pow(min_dist, 2))) / (3.*pow(max_dist - min_dist, 10));
	rb_coeffs[8][9] = (-3230 * sqrt(161)*pow(min_dist, 2)*(29 * max_dist + 26 * min_dist)) / (3.*pow(max_dist - min_dist, 10));
	rb_coeffs[8][10] = (17765 * sqrt(161)*pow(min_dist, 2)) / (3.*pow(max_dist - min_dist, 10));
	rb_coeffs[8][11] = 0;
	rb_coeffs[8][12] = 0;
	rb_coeffs[8][13] = 0;
	rb_coeffs[9][0] = -(sqrt(3)*pow(max_dist, 2)*pow(min_dist, 2)*(5 * pow(max_dist, 9) + 315 * pow(max_dist, 8)*min_dist + 5460 * pow(max_dist, 7)*pow(min_dist, 2) + 38220 * pow(max_dist, 6)*pow(min_dist, 3) + 126126 * pow(max_dist, 5)*pow(min_dist, 4) + 210210 * pow(max_dist, 4)*pow(min_dist, 5) + 180180 * pow(max_dist, 3)*pow(min_dist, 6) + 77220 * pow(max_dist, 2)*pow(min_dist, 7) + 15015 * max_dist*pow(min_dist, 8) + 1001 * pow(min_dist, 9))) / (4.*pow(max_dist - min_dist, 11));
	rb_coeffs[9][1] = (sqrt(3)*max_dist*pow(min_dist, 2)*(185 * pow(max_dist, 9) + 7035 * pow(max_dist, 8)*min_dist + 81900 * pow(max_dist, 7)*pow(min_dist, 2) + 405132 * pow(max_dist, 6)*pow(min_dist, 3) + 966966 * pow(max_dist, 5)*pow(min_dist, 4) + 1171170 * pow(max_dist, 4)*pow(min_dist, 5) + 720720 * pow(max_dist, 3)*pow(min_dist, 6) + 214500 * pow(max_dist, 2)*pow(min_dist, 7) + 27027 * max_dist*pow(min_dist, 8) + 1001 * pow(min_dist, 9))) / (2.*pow(max_dist - min_dist, 11));
	rb_coeffs[9][2] = -(sqrt(3)*pow(min_dist, 2)*(8885 * pow(max_dist, 9) + 227115 * pow(max_dist, 8)*min_dist + 1870596 * pow(max_dist, 7)*pow(min_dist, 2) + 6703788 * pow(max_dist, 6)*pow(min_dist, 3) + 11657646 * pow(max_dist, 5)*pow(min_dist, 4) + 10180170 * pow(max_dist, 4)*pow(min_dist, 5) + 4384380 * pow(max_dist, 3)*pow(min_dist, 6) + 859716 * pow(max_dist, 2)*pow(min_dist, 7) + 63063 * max_dist*pow(min_dist, 8) + 1001 * pow(min_dist, 9))) / (4.*pow(max_dist - min_dist, 11));
	rb_coeffs[9][3] = (6 * sqrt(3)*pow(min_dist, 2)*(4265 * pow(max_dist, 8) + 77196 * pow(max_dist, 7)*min_dist + 461188 * pow(max_dist, 6)*pow(min_dist, 2) + 1206296 * pow(max_dist, 5)*pow(min_dist, 3) + 1516515 * pow(max_dist, 4)*pow(min_dist, 4) + 930930 * pow(max_dist, 3)*pow(min_dist, 5) + 266266 * pow(max_dist, 2)*pow(min_dist, 6) + 30888 * max_dist*pow(min_dist, 7) + 1001 * pow(min_dist, 8))) / pow(max_dist - min_dist, 11);
	rb_coeffs[9][4] = (-51 * sqrt(3)*pow(min_dist, 2)*(3274 * pow(max_dist, 7) + 43022 * pow(max_dist, 6)*min_dist + 187824 * pow(max_dist, 5)*pow(min_dist, 2) + 355810 * pow(max_dist, 4)*pow(min_dist, 3) + 315315 * pow(max_dist, 3)*pow(min_dist, 4) + 129129 * pow(max_dist, 2)*pow(min_dist, 5) + 22022 * max_dist*pow(min_dist, 6) + 1144 * pow(min_dist, 7))) / pow(max_dist - min_dist, 11);
	rb_coeffs[9][5] = (714 * sqrt(3)*pow(min_dist, 2)*(942 * pow(max_dist, 6) + 9054 * pow(max_dist, 5)*min_dist + 28665 * pow(max_dist, 4)*pow(min_dist, 2) + 38350 * pow(max_dist, 3)*pow(min_dist, 3) + 22737 * pow(max_dist, 2)*pow(min_dist, 4) + 5577 * max_dist*pow(min_dist, 5) + 429 * pow(min_dist, 6))) / pow(max_dist - min_dist, 11);
	rb_coeffs[9][6] = (-6783 * sqrt(3)*pow(min_dist, 2)*(258 * pow(max_dist, 5) + 1800 * pow(max_dist, 4)*min_dist + 4030 * pow(max_dist, 3)*pow(min_dist, 2) + 3614 * pow(max_dist, 2)*pow(min_dist, 3) + 1287 * max_dist*pow(min_dist, 4) + 143 * pow(min_dist, 5))) / pow(max_dist - min_dist, 11);
	rb_coeffs[9][7] = (1938 * sqrt(3)*pow(min_dist, 2)*(1545 * pow(max_dist, 4) + 7630 * pow(max_dist, 3)*min_dist + 11466 * pow(max_dist, 2)*pow(min_dist, 2) + 6188 * max_dist*pow(min_dist, 3) + 1001 * pow(min_dist, 4))) / pow(max_dist - min_dist, 11);
	rb_coeffs[9][8] = (-969 * sqrt(3)*pow(min_dist, 2)*(6905 * pow(max_dist, 3) + 22911 * pow(max_dist, 2)*min_dist + 20748 * max_dist*pow(min_dist, 2) + 5096 * pow(min_dist, 3))) / (2.*pow(max_dist - min_dist, 11));
	rb_coeffs[9][9] = (3553 * sqrt(3)*pow(min_dist, 2)*(661 * pow(max_dist, 2) + 1323 * max_dist*min_dist + 546 * pow(min_dist, 2))) / pow(max_dist - min_dist, 11);
	rb_coeffs[9][10] = (-81719 * sqrt(3)*pow(min_dist, 2)*(23 * max_dist + 21 * min_dist)) / (2.*pow(max_dist - min_dist, 11));
	rb_coeffs[9][11] = (163438 * sqrt(3)*pow(min_dist, 2)) / pow(max_dist - min_dist, 11);
	rb_coeffs[9][12] = 0;
	rb_coeffs[9][13] = 0;
	rb_coeffs[10][0] = (3 * sqrt(0.5454545454545454)*pow(max_dist, 2)*pow(min_dist, 2)*(pow(max_dist, 10) + 75 * pow(max_dist, 9)*min_dist + 1575 * pow(max_dist, 8)*pow(min_dist, 2) + 13650 * pow(max_dist, 7)*pow(min_dist, 3) + 57330 * pow(max_dist, 6)*pow(min_dist, 4) + 126126 * pow(max_dist, 5)*pow(min_dist, 5) + 150150 * pow(max_dist, 4)*pow(min_dist, 6) + 96525 * pow(max_dist, 3)*pow(min_dist, 7) + 32175 * pow(max_dist, 2)*pow(min_dist, 8) + 5005 * max_dist*pow(min_dist, 9) + 273 * pow(min_dist, 10))) / pow(max_dist - min_dist, 12);
	rb_coeffs[10][1] = (-9 * sqrt(0.5454545454545454)*max_dist*pow(min_dist, 2)*(29 * pow(max_dist, 10) + 1325 * pow(max_dist, 9)*min_dist + 18900 * pow(max_dist, 8)*pow(min_dist, 2) + 117390 * pow(max_dist, 7)*pow(min_dist, 3) + 363090 * pow(max_dist, 6)*pow(min_dist, 4) + 594594 * pow(max_dist, 5)*pow(min_dist, 5) + 525525 * pow(max_dist, 4)*pow(min_dist, 6) + 246675 * pow(max_dist, 3)*pow(min_dist, 7) + 57915 * pow(max_dist, 2)*pow(min_dist, 8) + 5915 * max_dist*pow(min_dist, 9) + 182 * pow(min_dist, 10))) / pow(max_dist - min_dist, 12);
	rb_coeffs[10][2] = (9 * sqrt(0.5454545454545454)*pow(min_dist, 2)*(822 * pow(max_dist, 10) + 25525 * pow(max_dist, 9)*min_dist + 261135 * pow(max_dist, 8)*pow(min_dist, 2) + 1195740 * pow(max_dist, 7)*pow(min_dist, 3) + 2757300 * pow(max_dist, 6)*pow(min_dist, 4) + 3360357 * pow(max_dist, 5)*pow(min_dist, 5) + 2177175 * pow(max_dist, 4)*pow(min_dist, 6) + 725010 * pow(max_dist, 3)*pow(min_dist, 7) + 113490 * pow(max_dist, 2)*pow(min_dist, 8) + 6825 * max_dist*pow(min_dist, 9) + 91 * pow(min_dist, 10))) / pow(max_dist - min_dist, 12);
	rb_coeffs[10][3] = (-255 * sqrt(0.5454545454545454)*pow(min_dist, 2)*(397 * pow(max_dist, 9) + 8847 * pow(max_dist, 8)*min_dist + 66780 * pow(max_dist, 7)*pow(min_dist, 2) + 228228 * pow(max_dist, 6)*pow(min_dist, 3) + 392301 * pow(max_dist, 5)*pow(min_dist, 4) + 351351 * pow(max_dist, 4)*pow(min_dist, 5) + 162162 * pow(max_dist, 3)*pow(min_dist, 6) + 36270 * pow(max_dist, 2)*pow(min_dist, 7) + 3393 * max_dist*pow(min_dist, 8) + 91 * pow(min_dist, 9))) / pow(max_dist - min_dist, 12);
	rb_coeffs[10][4] = (2295 * sqrt(0.5454545454545454)*pow(min_dist, 2)*(345 * pow(max_dist, 8) + 5676 * pow(max_dist, 7)*min_dist + 32004 * pow(max_dist, 6)*pow(min_dist, 2) + 81627 * pow(max_dist, 5)*pow(min_dist, 3) + 103285 * pow(max_dist, 4)*pow(min_dist, 4) + 66066 * pow(max_dist, 3)*pow(min_dist, 5) + 20566 * pow(max_dist, 2)*pow(min_dist, 6) + 2769 * max_dist*pow(min_dist, 7) + 117 * pow(min_dist, 8))) / pow(max_dist - min_dist, 12);
	rb_coeffs[10][5] = (-8721 * sqrt(0.5454545454545454)*pow(min_dist, 2)*(444 * pow(max_dist, 7) + 5460 * pow(max_dist, 6)*min_dist + 22995 * pow(max_dist, 5)*pow(min_dist, 2) + 43225 * pow(max_dist, 4)*pow(min_dist, 3) + 39130 * pow(max_dist, 3)*pow(min_dist, 4) + 16926 * pow(max_dist, 2)*pow(min_dist, 5) + 3185 * max_dist*pow(min_dist, 6) + 195 * pow(min_dist, 7))) / pow(max_dist - min_dist, 12);
	rb_coeffs[10][6] = (20349 * sqrt(0.5454545454545454)*pow(min_dist, 2)*(612 * pow(max_dist, 6) + 5625 * pow(max_dist, 5)*min_dist + 17475 * pow(max_dist, 4)*pow(min_dist, 2) + 23530 * pow(max_dist, 3)*pow(min_dist, 3) + 14430 * pow(max_dist, 2)*pow(min_dist, 4) + 3783 * max_dist*pow(min_dist, 5) + 325 * pow(min_dist, 6))) / pow(max_dist - min_dist, 12);
	rb_coeffs[10][7] = (-8721 * sqrt(0.5454545454545454)*pow(min_dist, 2)*(3099 * pow(max_dist, 5) + 21025 * pow(max_dist, 4)*min_dist + 46830 * pow(max_dist, 3)*pow(min_dist, 2) + 42770 * pow(max_dist, 2)*pow(min_dist, 3) + 15925 * max_dist*pow(min_dist, 4) + 1911 * pow(min_dist, 5))) / pow(max_dist - min_dist, 12);
	rb_coeffs[10][8] = (43605 * sqrt(0.5454545454545454)*pow(min_dist, 2)*(913 * pow(max_dist, 4) + 4444 * pow(max_dist, 3)*min_dist + 6720 * pow(max_dist, 2)*pow(min_dist, 2) + 3731 * max_dist*pow(min_dist, 3) + 637 * pow(min_dist, 4))) / pow(max_dist - min_dist, 12);
	rb_coeffs[10][9] = (-111435 * sqrt(0.5454545454545454)*pow(min_dist, 2)*(352 * pow(max_dist, 3) + 1164 * pow(max_dist, 2)*min_dist + 1071 * max_dist*pow(min_dist, 2) + 273 * pow(min_dist, 3))) / pow(max_dist - min_dist, 12);
	rb_coeffs[10][10] = (334305 * sqrt(0.5454545454545454)*pow(min_dist, 2)*(74 * pow(max_dist, 2) + 149 * max_dist*min_dist + 63 * pow(min_dist, 2))) / pow(max_dist - min_dist, 12);
	rb_coeffs[10][11] = (-334305 * sqrt(0.5454545454545454)*pow(min_dist, 2)*(27 * max_dist + 25 * min_dist)) / pow(max_dist - min_dist, 12);
	rb_coeffs[10][12] = (1448655 * sqrt(0.5454545454545454)*pow(min_dist, 2)) / pow(max_dist - min_dist, 12);
	rb_coeffs[10][13] = 0;
	rb_coeffs[11][0] = -(sqrt(82.16666666666667)*pow(max_dist, 2)*pow(min_dist, 2)*(pow(max_dist, 11) + 88 * pow(max_dist, 10)*min_dist + 2200 * pow(max_dist, 9)*pow(min_dist, 2) + 23100 * pow(max_dist, 8)*pow(min_dist, 3) + 120120 * pow(max_dist, 7)*pow(min_dist, 4) + 336336 * pow(max_dist, 6)*pow(min_dist, 5) + 528528 * pow(max_dist, 5)*pow(min_dist, 6) + 471900 * pow(max_dist, 4)*pow(min_dist, 7) + 235950 * pow(max_dist, 3)*pow(min_dist, 8) + 62920 * pow(max_dist, 2)*pow(min_dist, 9) + 8008 * max_dist*pow(min_dist, 10) + 364 * pow(min_dist, 11))) / (4.*pow(max_dist - min_dist, 13));
	rb_coeffs[11][1] = (sqrt(82.16666666666667)*max_dist*pow(min_dist, 2)*(101 * pow(max_dist, 11) + 5456 * pow(max_dist, 10)*min_dist + 93500 * pow(max_dist, 9)*pow(min_dist, 2) + 711480 * pow(max_dist, 8)*pow(min_dist, 3) + 2762760 * pow(max_dist, 7)*pow(min_dist, 4) + 5861856 * pow(max_dist, 6)*pow(min_dist, 5) + 7002996 * pow(max_dist, 5)*pow(min_dist, 6) + 4719000 * pow(max_dist, 4)*pow(min_dist, 7) + 1746030 * pow(max_dist, 3)*pow(min_dist, 8) + 331760 * pow(max_dist, 2)*pow(min_dist, 9) + 28028 * max_dist*pow(min_dist, 10) + 728 * pow(min_dist, 11))) / (4.*pow(max_dist - min_dist, 13));
	rb_coeffs[11][2] = -(sqrt(82.16666666666667)*pow(min_dist, 2)*(1667 * pow(max_dist, 11) + 61754 * pow(max_dist, 10)*min_dist + 767360 * pow(max_dist, 9)*pow(min_dist, 2) + 4363590 * pow(max_dist, 8)*pow(min_dist, 3) + 12852840 * pow(max_dist, 7)*pow(min_dist, 4) + 20762742 * pow(max_dist, 6)*pow(min_dist, 5) + 18762744 * pow(max_dist, 5)*pow(min_dist, 6) + 9390810 * pow(max_dist, 4)*pow(min_dist, 7) + 2492490 * pow(max_dist, 3)*pow(min_dist, 8) + 318890 * pow(max_dist, 2)*pow(min_dist, 9) + 16016 * max_dist*pow(min_dist, 10) + 182 * pow(min_dist, 11))) / (2.*pow(max_dist - min_dist, 13));
	rb_coeffs[11][3] = (11 * sqrt(739.5)*pow(min_dist, 2)*(809 * pow(max_dist, 10) + 21740 * pow(max_dist, 9)*min_dist + 201990 * pow(max_dist, 8)*pow(min_dist, 2) + 871920 * pow(max_dist, 7)*pow(min_dist, 3) + 1957410 * pow(max_dist, 6)*pow(min_dist, 4) + 2395484 * pow(max_dist, 5)*pow(min_dist, 5) + 1611610 * pow(max_dist, 4)*pow(min_dist, 6) + 580840 * pow(max_dist, 3)*pow(min_dist, 7) + 104520 * pow(max_dist, 2)*pow(min_dist, 8) + 8060 * max_dist*pow(min_dist, 9) + 182 * pow(min_dist, 10))) / (2.*pow(max_dist - min_dist, 13));
	rb_coeffs[11][4] = (-1045 * sqrt(739.5)*pow(min_dist, 2)*(157 * pow(max_dist, 9) + 3156 * pow(max_dist, 8)*min_dist + 22272 * pow(max_dist, 7)*pow(min_dist, 2) + 73332 * pow(max_dist, 6)*pow(min_dist, 3) + 124852 * pow(max_dist, 5)*pow(min_dist, 4) + 113932 * pow(max_dist, 4)*pow(min_dist, 5) + 55328 * pow(max_dist, 3)*pow(min_dist, 6) + 13572 * pow(max_dist, 2)*pow(min_dist, 7) + 1482 * max_dist*pow(min_dist, 8) + 52 * pow(min_dist, 9))) / (4.*pow(max_dist - min_dist, 13));
	rb_coeffs[11][5] = (209 * sqrt(739.5)*pow(min_dist, 2)*(4569 * pow(max_dist, 8) + 69792 * pow(max_dist, 7)*min_dist + 375900 * pow(max_dist, 6)*pow(min_dist, 2) + 939400 * pow(max_dist, 5)*pow(min_dist, 3) + 1193920 * pow(max_dist, 4)*pow(min_dist, 4) + 787696 * pow(max_dist, 3)*pow(min_dist, 5) + 260988 * pow(max_dist, 2)*pow(min_dist, 6) + 39000 * max_dist*pow(min_dist, 7) + 1950 * pow(min_dist, 8))) / (4.*pow(max_dist - min_dist, 13));
	rb_coeffs[11][6] = (-1463 * sqrt(739.5)*pow(min_dist, 2)*(633 * pow(max_dist, 7) + 7383 * pow(max_dist, 6)*min_dist + 30200 * pow(max_dist, 5)*pow(min_dist, 2) + 56385 * pow(max_dist, 4)*pow(min_dist, 3) + 51870 * pow(max_dist, 3)*pow(min_dist, 4) + 23387 * pow(max_dist, 2)*pow(min_dist, 5) + 4732 * max_dist*pow(min_dist, 6) + 325 * pow(min_dist, 7))) / pow(max_dist - min_dist, 13);
	rb_coeffs[11][7] = (209 * sqrt(739.5)*pow(min_dist, 2)*(11814 * pow(max_dist, 6) + 104698 * pow(max_dist, 5)*min_dist + 320155 * pow(max_dist, 4)*pow(min_dist, 2) + 433020 * pow(max_dist, 3)*pow(min_dist, 3) + 272545 * pow(max_dist, 2)*pow(min_dist, 4) + 75166 * max_dist*pow(min_dist, 5) + 7007 * pow(min_dist, 6))) / pow(max_dist - min_dist, 13);
	rb_coeffs[11][8] = (-4807 * sqrt(739.5)*pow(min_dist, 2)*(3817 * pow(max_dist, 5) + 25300 * pow(max_dist, 4)*min_dist + 56080 * pow(max_dist, 3)*pow(min_dist, 2) + 51940 * pow(max_dist, 2)*pow(min_dist, 3) + 20020 * max_dist*pow(min_dist, 4) + 2548 * pow(min_dist, 5))) / (4.*pow(max_dist - min_dist, 13));
	rb_coeffs[11][9] = (24035 * sqrt(82.16666666666667)*pow(min_dist, 2)*(2959 * pow(max_dist, 4) + 14224 * pow(max_dist, 3)*min_dist + 21604 * pow(max_dist, 2)*pow(min_dist, 2) + 12264 * max_dist*pow(min_dist, 3) + 2184 * pow(min_dist, 4))) / (4.*pow(max_dist - min_dist, 13));
	rb_coeffs[11][10] = (-24035 * sqrt(82.16666666666667)*pow(min_dist, 2)*(1303 * pow(max_dist, 3) + 4294 * pow(max_dist, 2)*min_dist + 4000 * max_dist*pow(min_dist, 2) + 1050 * pow(min_dist, 3))) / (2.*pow(max_dist - min_dist, 13));
	rb_coeffs[11][11] = (28405 * sqrt(82.16666666666667)*pow(min_dist, 2)*(631 * pow(max_dist, 2) + 1276 * max_dist*min_dist + 550 * pow(min_dist, 2))) / (2.*pow(max_dist - min_dist, 13));
	rb_coeffs[11][12] = (-85215 * sqrt(739.5)*pow(min_dist, 2)*(47 * max_dist + 44 * min_dist)) / (4.*pow(max_dist - min_dist, 13));
	rb_coeffs[11][13] = (596505 * sqrt(739.5)*pow(min_dist, 2)) / (4.*pow(max_dist - min_dist, 13));
}

RadialBasis_Shapeev::RadialBasis_Shapeev(double _min_dist, double _max_dist, int _size)
	: AnyRadialBasis(_min_dist, _max_dist, _size)
{
	if (rb_size > ALLOCATED_DEGREE + 1)
		ERROR("RadialBasis error: allocated degree ecceded.");
	InitShapeevRB();
}

RadialBasis_Shapeev::RadialBasis_Shapeev(std::ifstream & ifs)
	: AnyRadialBasis(ifs)
{
	if (rb_size > ALLOCATED_DEGREE + 1)
		ERROR("RadialBasis error: allocated degree ecceded.");
	InitShapeevRB();
}

void RadialBasis_Shapeev::RB_Calc(double r)
{
#ifdef MLIP_DEBUG
	if (r < min_dist) {
		Warning("RadialBasis: r<min_dist. r = " + to_string(r) +
			", min_dist = " + to_string(min_dist) + '\n');
	}
	if (r > max_dist) {
		ERROR("RadialBasis: r>MaxDist !!!. r = " + to_string(r) +
			", min_dist = " + to_string(min_dist) + '\n');
	}
#endif

	for (int xi = 0; xi < rb_size; xi++) {
		rb_vals[xi] = 0;
		rb_ders[xi] = 0;
		for (int i = -2; i <= xi; i++) {
			rb_vals[xi] += scaling * rb_coeffs[xi][i + 2] * pow(r, i);
			rb_ders[xi] += scaling * rb_coeffs[xi][i + 2] * i * pow(r, i - 1);
		}
	}
}

void RadialBasis_Chebyshev::RB_Calc(double r)
{
#ifdef MLIP_DEBUG
	if (r < min_dist) {
		Warning("RadialBasis: r<min_dist. r = " + to_string(r) +
			", min_dist = " + to_string(min_dist) + '\n');
	}
	if (r > max_dist) {
		ERROR("RadialBasis: r>MaxDist !!!. r = " + to_string(r) +
			", min_dist = " + to_string(min_dist) + '\n');
	}
#endif

	double mult = 2.0 / (max_dist - min_dist);
	double ksi = (2 * r - (min_dist + max_dist)) / (max_dist - min_dist);

	rb_vals[0] = scaling * (1 * (r - max_dist)*(r - max_dist));
	rb_ders[0] = scaling * (0 * (r - max_dist)*(r - max_dist) + 2 * (r - max_dist));
	rb_vals[1] = scaling * (ksi*(r - max_dist)*(r - max_dist));
	rb_ders[1] = scaling * (mult * (r - max_dist)*(r - max_dist) + 2 * ksi*(r - max_dist));
	for (int i = 2; i < rb_size; i++) {
		rb_vals[i] = 2 * ksi*rb_vals[i - 1] - rb_vals[i - 2];
		rb_ders[i] = 2 * (mult * rb_vals[i - 1] + ksi * rb_ders[i - 1]) - rb_ders[i - 2];
	}
}

void RadialBasis_Chebyshev_repuls::RB_Calc(double r)
{
#ifdef MLIP_DEBUG
	if (r < min_dist) {
		Warning("RadialBasis: r<min_dist. r = " + to_string(r) +
			", min_dist = " + to_string(min_dist) + '\n');
	}
	if (r > max_dist) {
		ERROR("RadialBasis: r>MaxDist !!!. r = " + to_string(r) +
			", min_dist = " + to_string(min_dist) + '\n');
	}
#endif

	if (r < min_dist)
		r = min_dist;

	double mult = 2.0 / (max_dist - min_dist);
	double ksi = (2 * r - (min_dist + max_dist)) / (max_dist - min_dist);

	rb_vals[0] = scaling * (1 * (r - max_dist)*(r - max_dist));
	rb_ders[0] = scaling * (0 * (r - max_dist)*(r - max_dist) + 2 * (r - max_dist));
	rb_vals[1] = scaling * (ksi*(r - max_dist)*(r - max_dist));
	rb_ders[1] = scaling * (mult * (r - max_dist)*(r - max_dist) + 2 * ksi*(r - max_dist));
	for (int i = 2; i < rb_size; i++) {
		rb_vals[i] = 2 * ksi*rb_vals[i - 1] - rb_vals[i - 2];
		rb_ders[i] = 2 * (mult * rb_vals[i - 1] + ksi * rb_ders[i - 1]) - rb_ders[i - 2];
	}
	if (r == min_dist)
		for (int i = 0; i < rb_size; i++)
			rb_ders[i] = 0.0;
}

void RadialBasis_Taylor::RB_Calc(double r)
{
#ifdef MLIP_DEBUG
	if (r < min_dist) {
		Warning("RadialBasis: r<min_dist. r = " + to_string(r) +
			", min_dist = " + to_string(min_dist) + '\n');
	}
	if (r > max_dist) {
		ERROR("RadialBasis: r>MaxDist !!!. r = " + to_string(r) +
			", min_dist = " + to_string(min_dist) + '\n');
	}
#endif

	rb_vals[0] = scaling * 1;
	rb_ders[0] = scaling * 0;
	for (int i = 1; i < rb_size; i++)
	{
		rb_ders[i] = scaling * i * rb_vals[i - 1];
		rb_vals[i] = scaling * r * rb_vals[i - 1];
	}
}

const double RadialBasis_vdW::RECURSIVE_COEFF[22][3] = { {0,0,0}, {8.87248718788328392,-0.542267396590896327,0},{6.79251268298387606,-0.509113772340205375,-0.112707967768682771},{5.97200240050575452,-0.490871142415655530,-0.147220924961263526},{5.53059832007591211,-0.480355753378264948,-0.167448023784336122},{5.25531942342471711,-0.474098576668924181,-0.180812263362180705},{5.06759630263062248,-0.470320600452479026,-0.190283390871098224},{4.93157096163622518,-0.468051693398248887,-0.197332214383551713},{4.82855147527615044,-0.466734005666727749,-0.202775141588597198},{4.74785129937755951,-0.466032246548252509,-0.207101447529418508},{4.68292754573870735,-0.465737360308779940,-0.210621592157087861},{4.62955581831805592,-0.465714968307180773,-0.213541633995589659},{4.58489017048410400,-0.465876479156477367,-0.216003443795458052},{4.54694521651720235,-0.466162247521016685,-0.218107732751734198},{4.51429486932813924,-0.466531414516049284,-0.219927875173733517},{4.48588903210495412,-0.466955593424251038,-0.221518538098693055},{4.46093775571628900,-0.467414841652198140,-0.222921252140461663},{4.43883562185532690,-0.467895031561814091,-0.224168113244483244},{4.41911098018499224,-0.468386099282248837,-0.225284305432789144},{4.40139102397022533,-0.468880857197526929,-0.226289858861643283},{4.38537723625471118,-0.469374175729228661,-0.227200899568782515},{4.37082779110205071,-0.469862411506819840,-0.228030553844448804} };
const double RadialBasis_vdW::ZEROTH_POLY = 102.295067549833082;

void RadialBasis_vdW::RB_Calc(double r) {

	
	const double r_sq = r*r;

	const double mindist_sq = min_dist * min_dist;
	const double mindist_sq_plus_eps = mindist_sq * 1.02;

	if (r_sq <= mindist_sq_plus_eps)
	{
		FillWithZero(rb_vals);
		rb_vals[0] = 2.5;
		FillWithZero(rb_ders);
		return;
	}
	if (r >= max_dist)
	{
		FillWithZero(rb_vals);
		FillWithZero(rb_ders);
		return;
	}


	const double x_inv_sq = r_sq / (min_dist*min_dist);
	const double x_sq = 1.0 / x_inv_sq;

	// cout << "x_sq:" << x_sq << endl;

	const double x_sq_minus_one_inv = 1.0 / (x_sq - 1.0);

	// if (x_inv_sq < 1.01)
	// 	cout << "x_sq_minus_one_inv:" << x_sq_minus_one_inv << endl;

	const double x_sq_minus_one_inv_sq = x_sq_minus_one_inv * x_sq_minus_one_inv;
	const double my_exp = exp(x_sq_minus_one_inv);

	// if (x_inv_sq < 1.01)
	// 	cout << "my_exp" << my_exp << endl;

	const double x_sq_2 = x_sq * x_sq;
	const double x_sq_3 = x_sq_2 * x_sq;
	const double mult = x_sq_3 * my_exp;
	const double drinv_dr = -1.0 * x_sq / r_sq;

	double e_times_my_exp = 2.71828182845904524 * my_exp;
	rb_vals[0] = 2.5 * std::pow(1 - e_times_my_exp, 3);

	// if (x_inv_sq < 1.01)
	// 	cout << "rb_vals[0]:" << rb_vals[0] << endl;

	rb_ders[0] = 7.5 * std::pow(1 - e_times_my_exp, 2) * e_times_my_exp * x_sq_minus_one_inv_sq * drinv_dr;

	// if (x_inv_sq < 1.01)
	// 	cout << "rb_ders[0]:" << rb_ders[0] << endl;

	if (rb_size == 1) return;
	rb_vals[1] = ZEROTH_POLY * mult;
	rb_ders[1] = ZEROTH_POLY * my_exp * x_sq_2 * x_sq_minus_one_inv_sq * ((3 * x_sq - 7) * x_sq + 3) * drinv_dr;
	double prev = 0;
	double prev_der = 0;
	for (int i = 1; i < rb_size - 1; i++) {
		rb_vals[i + 1] = RECURSIVE_COEFF[i][0] * (
			(x_sq + RECURSIVE_COEFF[i][1]) * rb_vals[i] + RECURSIVE_COEFF[i][2] * prev);
		rb_ders[i + 1] = RECURSIVE_COEFF[i][0] * (
			drinv_dr * rb_vals[i] +
			(x_sq + RECURSIVE_COEFF[i][1]) * rb_ders[i] + RECURSIVE_COEFF[i][2] * prev_der);
		prev = rb_vals[i];
		prev_der = rb_ders[i];
	}

	for (int i = 0; i < rb_size; i++)
	{
		rb_ders[i] *= 2*r;
	}
}
