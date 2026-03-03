#pragma once
#include <iostream>
#include <array>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <sstream>
#include <random>
#include <map>
#include <unordered_map>
#include <utility>
#include <algorithm>
#include <chrono>
#include <iomanip>

/* DEFINE DATA TYPES FOR CONVENIENCE -------------------------------------- */
using UShort = unsigned short;

class Timer
{
	private:
		std::chrono::time_point<std::chrono::steady_clock> m_beg{std::chrono::steady_clock::now()};
	public:
		void reset() { m_beg = std::chrono::steady_clock::now(); }
		float elapsed() const
		{
			return std::chrono::duration_cast<std::chrono::duration<float,std::ratio<1>>> (std::chrono::steady_clock::now()-m_beg).count();
		}
};

/* CUSTOM TYPES ------------------------------------------------------------- */
struct ATOM
{
	unsigned int					total;
	std::vector<uint>				ID, NET;
	std::vector<float>				Z;
	std::vector<double>				X, Y;
	std::vector<bool>				PHANTOM;
	std::vector<std::array<uint,6>>	BOND;
	std::vector<std::array<bool,4>>	BRTL;	
}
;

struct BOND_NET
{
	uint				total;
	std::vector<UShort> PRD;
	std::vector<uint> 	ID, A, B;
	std::vector<double> L0;
};

struct ANGLE_NET
{
	uint				total;
	std::vector<UShort> PRD1, PRD2;
	std::vector<uint> 	ID, A, B, C;
	std::vector<double> T0;
};

struct POINT
{
	uint	total;
	uint	*ID;
	double	*X, *Y, *FX, *FY, *VX, *VY;

	void allocate_memory()
	{
		ID = new uint[total];
		X = new double[total];	Y = new double[total];
		VX = new double[total];	VY = new double[total];
		FX = new double[total];	FY = new double[total];
	}

	void free_memory()
	{
		delete[] ID;
		delete[] X;		delete[] Y;
		delete[] VX;	delete[] VY;
		delete[] FX;	delete[] FY;
	}
};

struct BOND
{
	uint	total;
	UShort	*PRD;
	uint	*ID, *A, *B;
	double	*L0;

	void allocate_memory()
	{
		ID = new uint[total];
		A = new uint[total]; B = new uint[total];
		L0 = new double[total];
		PRD = new UShort[total];
	}
};

struct ANGLE
{
	uint	total;
	UShort	*PRD1, *PRD2;
	uint	*ID, *A, *B, *C;
	double	*T0;

	void allocate_memory()
	{
		ID = new uint[total];
		A = new uint[total]; B = new uint[total]; C = new uint[total];
		T0 = new double[total];
		PRD1 = new UShort[total]; PRD2 = new UShort[total];
	}
};

struct SHIFT_BOND
{
	double	*ax, *ay, *bx, *by;

	void allocate_memory( const int& size )
	{
		ax = new double[size];	ay = new double[size];
		bx = new double[size];	by = new double[size];
	}
};

struct SHIFT_ANGLE
{
	double	*ax, *ay, *bx1, *by1, *bx2, *by2, *cx, *cy;

	void allocate_memory( const int& size )
	{
		ax = new double[size];	ay = new double[size];
		bx1 = new double[size];	by1 = new double[size];
		bx2 = new double[size];	by2 = new double[size];
		cx = new double[size];	cy = new double[size];
	}
};

/* RESERVE MEMORY ---------------------------------------------------------- */
void reserve_memory( ATOM& atoms, BOND_NET& bonds, ANGLE_NET& angles, const int& total_atoms )
{
	int atoms_size { total_atoms*5 }, angles_size { total_atoms*6 }, bonds_size { total_atoms*6 };

	atoms.ID.reserve( atoms_size );		angles.ID.reserve( angles_size );	bonds.ID.reserve( bonds_size );	
	atoms.NET.reserve( atoms_size );	angles.PRD1.reserve( angles_size );	bonds.PRD.reserve( bonds_size );
	atoms.X.reserve( atoms_size );		angles.PRD2.reserve( angles_size );	bonds.A.reserve( bonds_size );
	atoms.Y.reserve( atoms_size );		angles.A.reserve( angles_size );	bonds.B.reserve( bonds_size );
	atoms.Z.reserve( atoms_size );		angles.B.reserve( angles_size );	bonds.L0.reserve( bonds_size );
	atoms.BOND.reserve( atoms_size );	angles.C.reserve( angles_size );
	atoms.BRTL.reserve( atoms_size );	angles.T0.reserve( angles_size );
	atoms.PHANTOM.reserve( atoms_size );
}

/* LABEL CONNECTED CLUSTERS OF ATOMS ---------------------------------------- */
void set_network ( const uint& id, ATOM& atoms )
{
	uint temp;
	for ( int j = 0; j < 3; j++ )
	{
		temp = atoms.BOND[id][j];
		if ( temp && !atoms.NET[temp-1] )
		{
			atoms.NET[temp-1] = atoms.NET[id];
			set_network( temp-1, atoms );
		}
	}
	for ( int j = 3; j < 6; j++ )
	{
		temp = atoms.BOND[id][j];
		if ( temp )
		{
			if ( atoms.BOND[temp-1][j-3] && !atoms.NET[temp-1] )
			{
				atoms.NET[temp-1] = atoms.NET[id];
				set_network ( temp-1, atoms);
			}
		}
	}
}

/* CONNECTIVITY OF EACH INDIVIDUAL ATOM ------------------------------------------- */
void localZ ( ATOM& atoms )
{
	for ( int i = 0; i < atoms.total; i++ )
	{
		atoms.Z[i] = 0;
		for ( int j = 0; j < 6; j++ )	if ( atoms.BOND[i][j] )	atoms.Z[i]++;
	}
}

/* AVERAGE CONNECTIVITY OF NETWORK ------------------------------------------- */
float calcZ ( const ATOM& atoms )
{
	float	avgZ {}, count {};

	for ( int i = 0; i < atoms.total; i++ )
	{
		if ( !atoms.PHANTOM[i] )
		{
			avgZ += atoms.Z[i];
			count += 1.0f;
		}
	}

	return avgZ / count;
}

/* REMOVE ATOMS WITH Z=1 (RECURSIVE) ----------------------------------------------- */
void remove_dangling_ends ( uint id, ATOM& atoms )
{
	uint temp;

	for(int j=0; j<6; j++)
	{
		temp = atoms.BOND[id][j];
		if ( temp )
		{
			atoms.BOND[id][j] = 0;
			atoms.Z[id] -= 1;
			atoms.ID[id] = 0;

			atoms.BOND[temp-1][(j+3)%6] = 0;
			atoms.Z[temp-1] -= 1;
			if ( atoms.Z[temp-1] == 1 )	remove_dangling_ends ( temp-1, atoms );
			break;
		}
	}
}

/* SHORTEST CARTESIAN DISTANCE CHECK ----------------------------------------------------------- */
bool minDist ( double& dist1sq, double& ax, double& ay, double& bx, double& by )
{
	double dist2sq = (ax-bx)*(ax-bx) + (ay-by)*(ay-by);
	if ( dist2sq < dist1sq )
	{
		dist1sq = dist2sq;
		return true;
	}
	else return false;
}

/* LABEL EACH BOND ACCORDING TO BOUNDARY ------------------------------------------------- */
UShort periodic ( double& AX, double& AY, double& BX, double& BY, const std::array<bool,4>& brtlA, const std::array<bool,4>& brtlB, UShort& W )
{
	const double Wsqrt3by2 { W*0.866025403784439 };
	UShort	prd {};
	double	ax, ay, bx, by, dist1sq;

	ax = AX; ay = AY; bx = BX; by = BY;	// reset
	dist1sq = (ax-bx)*(ax-bx) + (ay-by)*(ay-by);

	if ( brtlA[0] && brtlB[2] )		// UP (a)
	{
		ay += Wsqrt3by2;
		if (minDist( dist1sq, ax, ay, bx, by )) prd=1;
	}
	else if ( brtlA[2] && brtlB[0] )	// UP (b)
	{
		by += Wsqrt3by2;
		if (minDist( dist1sq, ax, ay, bx, by )) prd=2;
	}

	ax = AX; ay = AY; bx = BX; by = BY;	// reset

	if ( brtlA[3] && brtlB[1] )		// RIGHT (a)
	{
		ax += W;
		if (minDist( dist1sq, ax, ay, bx, by )) prd=3;
	}
	else if ( brtlA[1] && brtlB[3] )	// RIGHT (b)
	{
		bx += W;
		if (minDist( dist1sq, ax, ay, bx, by )) prd=4;
	}

	ax = AX; ay = AY; bx = BX; by = BY;	// reset

	if ( brtlA[3] && brtlA[0] && brtlB[1] && brtlB[2] )	// UP (a)-RIGHT (a)
	{
		ax += W;
		ay += Wsqrt3by2;
		if (minDist( dist1sq, ax, ay, bx, by )) prd=5;
	}
	else if ( brtlA[1] && brtlA[2] && brtlB[3] && brtlB[0] )	// UP (b)-RIGHT (b)
	{
		bx += W;
		by += Wsqrt3by2;
		if (minDist( dist1sq, ax, ay, bx, by )) prd=6;
	}

	return prd;
}

/* CREATE TEMPORARY MIRROR POINTS (FOR ANGLES) ------------------------------------------------- */
template <typename T1>
void shiftPointsA ( const T1& angles, SHIFT_ANGLE& shiftA, const UShort& W, const double& strain, const double& deform )
{
	constexpr double sqrt3by2 { 0.866025403784439 };
	const double	 upX {W*sqrt3by2*(1.0+deform)*strain}, upY {W*(1.0+deform)*sqrt3by2}, upXplusW {W*(1.0+sqrt3by2*(1.0+deform)*strain)};

	for ( int i = 0; i < angles.total; ++i )
	{
		shiftA.ax[i] = 0.0;		shiftA.ay[i] = 0.0;
		shiftA.bx1[i] = 0.0;	shiftA.by1[i] = 0.0;
		shiftA.bx2[i] = 0.0;	shiftA.by2[i] = 0.0;
		shiftA.cx[i] = 0.0;		shiftA.cy[i] = 0.0;

		switch ( angles.PRD1[i] )
		{
			case 0:
				break;
			case 1:		// Up (a)
				shiftA.ax[i] = upX;
				shiftA.ay[i] = upY;
				break;
			case 2:		// Up (b)
				shiftA.bx1[i] = upX;
				shiftA.by1[i] = upY;
				break;
			case 3:		// Right (a)
				shiftA.ax[i] = W;
				break;
			case 4:		// Right (b)
				shiftA.bx1[i] = W;
				break;
			case 5:		// Up (a) & Right (a)
				shiftA.ax[i] = upXplusW;
				shiftA.ay[i] = upY;
				break;
			case 6:		// Up (b) & Right (b)
				shiftA.bx1[i] = upXplusW;
				shiftA.by1[i] = upY;
				break;
		}
		switch ( angles.PRD2[i] )
		{
			case 0:
				break;
			case 1:		// Up (a)
				shiftA.bx2[i] = upX;
				shiftA.by2[i] = upY;
				break;
			case 2:		// Up (b)
				shiftA.cx[i] = upX;
				shiftA.cy[i] = upY;
				break;
			case 3:		// Right (a)
				shiftA.bx2[i] = W;
				break;
			case 4:		// Right (b)
				shiftA.cx[i] = W;
				break;
			case 5:		// Up (a) & Right (a)
				shiftA.bx2[i] = upXplusW;
				shiftA.by2[i] = upY;
				break;
			case 6:		// Up (b) & Right (b)
				shiftA.cx[i] = upXplusW;
				shiftA.cy[i] = upY;
				break;
		}
	}
}

/* CREATE TEMPORARY MIRROR POINTS (FOR BONDS) ------------------------------------------------- */
template <typename T2>
void shiftPointsB ( const T2& bonds, SHIFT_BOND& shiftB, const UShort& W, const double& strain, const double& deform )
{
	constexpr double sqrt3by2 { 0.866025403784439 };
	const double	 upX {W*sqrt3by2*(1.0+deform)*strain}, upY {W*(1.0+deform)*sqrt3by2}, upXplusW {W*(1.0+sqrt3by2*(1.0+deform)*strain)};

	for ( int i = 0; i < bonds.total; ++i )
	{
		shiftB.ax[i] = 0.0;		shiftB.ay[i] = 0.0;
		shiftB.bx[i] = 0.0;		shiftB.by[i] = 0.0;

		switch ( bonds.PRD[i] )
		{
			case 0:
				break;
			case 1:		// Up (a)
				shiftB.ax[i] = upX;
				shiftB.ay[i] = upY;
				break;
			case 2:		// Up (b)
				shiftB.bx[i] = upX;
				shiftB.by[i] = upY;
				break;
			case 3:		// Right (a)
				shiftB.ax[i] = W;
				break;
			case 4:		// Right (b)
				shiftB.bx[i] = W;
				break;
			case 5:		// Up (a) & Right (a)
				shiftB.ax[i] = upXplusW;
				shiftB.ay[i] = upY;
				break;
			case 6:		// Up (b) & Right (b)
				shiftB.bx[i] = upXplusW;
				shiftB.by[i] = upY;
				break;
		}
	}
}

/* UPDATE DISTANCE ACCORDING TO PERIODICITY ------------------------------------------------- */
void nDist ( const ATOM& atoms, BOND_NET& bonds, SHIFT_BOND& shiftB, const UShort& W )
{
	double	dx, dy;
	shiftPointsB<BOND_NET>( bonds, shiftB, W, 0.0, 0.0 );

	for ( int i = 0; i < bonds.total; i++ )
	{
		dx = atoms.X[ bonds.A[i]-1 ] + shiftB.ax[i] - atoms.X[ bonds.B[i]-1 ] - shiftB.bx[i];
		dy = atoms.Y[ bonds.A[i]-1 ] + shiftB.ay[i] - atoms.Y[ bonds.B[i]-1 ] - shiftB.by[i];
		bonds.L0[i] = std::sqrt( dx*dx + dy*dy );
	}
}

/* UPDATE ANGLE ACCORDING TO PERIODITICITY ------------------------------------------------- */
void nDegrees ( const ATOM& atoms, ANGLE_NET& angles, SHIFT_ANGLE& shiftA, const UShort& W )
{
	double	theta, bax, bay, bcx, bcy;

	shiftPointsA<ANGLE_NET>( angles, shiftA, W, 0.0, 0.0 );

	for ( int i = 0; i < angles.total; i++ )
	{
		bax = atoms.X[angles.A[i]-1] + shiftA.ax[i] - atoms.X[angles.B[i]-1] - shiftA.bx1[i];
		bay = atoms.Y[angles.A[i]-1] + shiftA.ay[i] - atoms.Y[angles.B[i]-1] - shiftA.by1[i];
		bcx = atoms.X[angles.C[i]-1] + shiftA.cx[i] - atoms.X[angles.B[i]-1] - shiftA.bx2[i];
		bcy = atoms.Y[angles.C[i]-1] + shiftA.cy[i] - atoms.Y[angles.B[i]-1] - shiftA.by2[i];

		theta = (bax*bcx + bay*bcy) / std::sqrt( (bax*bax + bay*bay)*(bcx*bcx + bcy*bcy) );

		if (std::abs(theta) <= 1.0)	angles.T0[i] = std::acos(theta);	// radians
		else if (std::abs(theta) < 1.00001)
		{
			if (theta > 1.0)		theta = 1.0;
			else if (theta < -1.0)	theta = -1.0;
			angles.T0[i] = std::acos(theta);
		}
		else
		{
			std::cout << "ANGLE:\t" << angles.ID[i] << '\t' << theta << '\n';
			std::exit(1);
		}
	}
}