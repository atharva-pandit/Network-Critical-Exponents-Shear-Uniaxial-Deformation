#include "setup.hpp"
#include "read_write.hpp"
#include "minimise.hpp"

int main( int argc, char* argv[] )
{
	/* INPUT ----------------------------------------------------------------------- */
	constexpr UShort	steps {100};
	constexpr double	timestep {0.05}, firstStrain {1e-3}, finalStrain {10.0};

	/* VARIABLES ------------------------------------------------------------------ */
	float				Z;
	double				energy, stepsize, strain {}, deform, def {}, dirxn;
	UShort				W, K, V;
	POINT				points;
	BOND				bonds;
	ANGLE				angles;
	std::stringstream	arg1 {argv[1]}, arg2 {argv[2]}, arg3 {argv[3]}, arg4 {argv[4]}, arg5 {argv[5]}, arg6 {argv[6]};
	Timer				Time;

	/* DEFINE KAPPA, Z & EXTENSION USING TERMINAL INPUT ---------------------------------------- */
	arg1 >> W;	arg2 >> Z;	arg3 >> K;	arg4 >> deform;	arg5 >> V;	arg6 >> dirxn;
	const double kappa = std::pow( 10.0, -K );

	/* READ FILENAMES & WRITE ENERGY -------------------------------------------------------------------------------- */
	const std::string	pntfile {"data/p"+mkStr(W,0)+"_Z"+mkStr(Z,2)+"_V"+mkStr(V,0)+".dat"},
						bndfile {"data/b"+mkStr(W,0)+"_Z"+mkStr(Z,2)+"_V"+mkStr(V,0)+".dat"},
						angfile {"data/a"+mkStr(W,0)+"_Z"+mkStr(Z,2)+"_V"+mkStr(V,0)+".dat"},
						enefile {"stress/W"+mkStr(W,0)+"_Z"+mkStr(Z,2)+"_K"+mkStr(K,0)+"_E"+mkStr(deform,2)+"_V"+mkStr(V,0)+"_r"+mkStr(dirxn,0)+".dat"};

	/* READ LATTICE DATA & SET UP MIRROR POINTS ---------------------------------------------------------- */
	readPoint ( pntfile, points );
	readBond ( bndfile, bonds );
	readAngle ( angfile, angles );
	SHIFT_ANGLE	shiftA;		shiftA.allocate_memory( angles.total );
	SHIFT_BOND	shiftB;		shiftB.allocate_memory( bonds.total );

	/* LOGARITHMIC STRAIN ---------------------------------------------------------- */
	std::array<double, steps+1>	strain_array {};
	stepsize = std::log10(finalStrain / firstStrain) / (steps-1);
	for (int i=0; i<steps; ++i) strain_array[i+1] = std::pow(10.0, std::log10(firstStrain) + i*stepsize);

	/* PRINT INFO & CREATE/CLEAR STRESS-STRAIN FILE ---------------------------- */
	std::cout << std::fixed << std::setprecision(2) << "RUN\tW: " << W << "\tZ: " << Z << "\tV: " << V << "\tK: " << K << "\tE: " << deform << "\tR: " << dirxn << '\n';
	clearfile ( enefile );

	/* COMPRESS(-)/EXTEND(+) & RELAX BEFORE SHEARING -------------------------------------------------- */
 	const double increment = deform >= 0.0 ? 0.001 : -0.001;
	strain = 0.0;

	while ( (std::abs(deform) - std::abs(def)) > 1e-6 )
	{
		def += increment;
		for (int i=0; i<points.total; ++i)	points.Y[i] += points.Y[i] * increment;
		energy = FIRE ( points, bonds, angles, timestep, W, kappa, shiftA, shiftB, strain, def );
	}

	writeEnergy<double>( strain, energy, enefile );
	//saveState( points, W, Z, deform, K, V, dirxn, 0 );

	/* SHEAR & MINIMIZE ------------------------------------------------------------ */
 	for (int tstep=0; tstep<steps; tstep++)
	{
		// Apply incremental shear
		for (int i=0; i<points.total; ++i) points.X[i] += points.Y[i] * (strain_array[tstep+1] - strain_array[tstep]) * dirxn;

		// Minimise Energy
		strain = strain_array[tstep+1] * dirxn;
		energy = FIRE ( points, bonds, angles, timestep, W, kappa, shiftA, shiftB, strain, deform );
		writeEnergy<double>( strain, energy, enefile );
		//saveState( points, W, Z, deform, K, V, dirxn, tstep+1 );
	}

	std::cout << std::fixed << std::setprecision(2) <<
	"DONE\tW: " << W << "\tZ: " << Z << "\tV: " << V <<"\tK: " << K << "\tE: " << deform << "\tR: " << dirxn << "\tTime: " << Time.elapsed() << " sec.\n";

	return 0;
}