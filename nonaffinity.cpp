#include "setup.hpp"
#include "read_write.hpp"

/* CALCULATE NON-AFFINITY ------------------------------------------------------- */
double nonaffinity ( const POINT& points, const POINT& pointsOLD, const double& perturbation )
{
	double			dx, dy, NA {};

	for ( int i = 0; i < points.total; ++i )
	{
		dx = points.X[i] - pointsOLD.X[i];
		dy = points.Y[i] - pointsOLD.Y[i];
		
		NA += dx*dx + dy*dy;
	}

	NA /= ( points.total * perturbation * perturbation );

	return NA;
}

/* MAIN ------------------------------------------------------------------------- */
int main( int argc, char* argv[] )
{
	/* INPUT ----------------------------------------------------------------------- */
	constexpr UShort	steps {100};
	constexpr double	timestep {0.05}, firstStrain {1e-3}, finalStrain {10.0};

	/* VARIABLES ------------------------------------------------------------------ */
	float				Z;
	double				stepsize, perturbation, nonaff, deform, dirxn;
	UShort				W, K, V;
	POINT				pointsOLD, pointsNEW;
	std::string 		pntfileOLD, pntfileNEW;
	std::stringstream	arg1 {argv[1]}, arg2 {argv[2]}, arg3 {argv[3]}, arg4 {argv[4]}, arg5 {argv[5]}, arg6 {argv[6]};
	Timer				Time;

	/* DEFINE KAPPA, Z & EXTENSION USING TERMINAL INPUT ---------------------------------------- */
	arg1 >> W;	arg2 >> Z;	arg3 >> K;	arg4 >> deform;	arg5 >> V;	arg6 >> dirxn;

	/* IN & OUT FILENAMES -------------------------------------------------------------------------------- */
	const std::string	NAfile	{"nonaffine/W"+mkStr(W,0)+"_Z"+mkStr(Z,2)+"_K"+mkStr(K,0)+"_E"+mkStr(deform,2)+"_V"+mkStr(V,0)+"_r"+mkStr(dirxn,0)+".dat"},
						//mainfile {"save/p"+mkStr(W,0)+"_Z"+mkStr(Z,2)+"_K"+mkStr(K,0)+"_E"+mkStr(deform,2)+"_V"+mkStr(V,0)+"_r"+mkStr(dirxn,0)};
						mainfile {"/media/panditat/Crucial X8/save_big/p"+mkStr(W,0)+"_Z"+mkStr(Z,2)+"_K"+mkStr(K,0)+"_E"+mkStr(deform,2)+"_V"+mkStr(V,0)+"_r"+mkStr(dirxn,0)};

	/* LOGARITHMIC STRAIN ---------------------------------------------------------- */
	std::array<double, steps+1>	strain_array {};
	stepsize = std::log10(finalStrain / firstStrain) / (steps-1);
	for (int i=0; i<steps; ++i) strain_array[i+1] = std::pow(10.0, std::log10(firstStrain) + i*stepsize);

	/* PRINT INFO & CREATE/CLEAR STRESS-STRAIN FILE ---------------------------- */
	std::cout << std::fixed << std::setprecision(2) << "RUN\tW: " << W << "\tZ: " << Z << "\tV: " << V << "\tK: " << K << "\tE: " << deform << "\tR: " << dirxn << '\n';
	clearfile ( NAfile );

	/* SHEAR ------------------------------------------------------------------------------------------------ */
 	for (int tstep=0; tstep<steps; tstep++)
	{
		// Update points ( old & new )
		pntfileOLD = mainfile + "_t"+mkStr(tstep,0)+".dat";
		pntfileNEW = mainfile + "_t"+mkStr(tstep+1,0)+".dat";
		readPoint ( pntfileOLD, pointsOLD );
		readPoint ( pntfileNEW, pointsNEW );

		// Incremental shear
		perturbation = strain_array[tstep+1] - strain_array[tstep];
		for (int i=0; i<pointsOLD.total; ++i) pointsOLD.X[i] += pointsOLD.Y[i] * perturbation * dirxn;

		// Calculate rearrangement
		nonaff = nonaffinity( pointsNEW, pointsOLD, perturbation );
		writeEnergy<double>( strain_array[tstep+1], nonaff, NAfile );

		// Free memory
		pointsOLD.free_memory();
		pointsNEW.free_memory();
	}

	std::cout << std::fixed << std::setprecision(2) <<
	"DONE\tW: " << W << "\tZ: " << Z << "\tV: " << V <<"\tK: " << K << "\tE: " << deform << "\tR: " << dirxn << "\tTime: " << Time.elapsed() << " sec.\n";

	return 0;
}