/* CALCULATE BOND LENGTH & UPDATE SPRING FORCE ----------------------------------------------------------- */
double updateForce_stretch_spring ( POINT& points, const BOND& bonds, const SHIFT_BOND& shiftB )
{
	uint	A, B;
	double	dx, dy,	forceX, forceY, bond_len, extensn, extensn_norm, energy {};

	for ( int i = 0; i < bonds.total; ++i )
	{
		A = bonds.A[i];	B = bonds.B[i];

		dx = points.X[B] + shiftB.bx[i] - points.X[A] - shiftB.ax[i];
		dy = points.Y[B] + shiftB.by[i] - points.Y[A] - shiftB.ay[i];

		bond_len = std::sqrt( dx*dx + dy*dy );
		extensn = bond_len - bonds.L0[i];
		extensn_norm = extensn / bond_len;

		forceX = extensn_norm * dx;
		forceY = extensn_norm * dy;

 		points.FX[A] += forceX;
		points.FY[A] += forceY;
		points.FX[B] -= forceX;
		points.FY[B] -= forceY;

		energy += extensn * extensn;
	}
	return	energy;
}

/* UPDATE FORCE ON EACH ATOM ACC TO CURRENT BOND LENGHTS & ANGLES ---------------------------------------------------------------- */
double updateForce_spring ( POINT& points, const BOND& bonds, const SHIFT_BOND& shiftB )
{
	double	energy_stretch, energy_bend;
	for ( int i = 0; i < points.total; ++i )	{ points.FX[i] = 0; points.FY[i] = 0; }

	energy_stretch = updateForce_stretch_spring( points, bonds, shiftB );

	return 0.5 * energy_stretch;
}

/* FIND MAXIMUM FORCE AT ANY NODE ------------------------------------------------------------------------- */
bool maxForce_spring ( const POINT& points )
{
	for (int i = 0; i<points.total; i++) if (std::sqrt(points.FX[i]*points.FX[i] + points.FY[i]*points.FY[i]) > 1e-12) return false;
	return true;
}

/* CALCULATE VIRIAL STRESS ------------------------------------------------------------------------------------------------- */
double virialStress_spring( const POINT& points, const BOND& bonds, const SHIFT_BOND& shiftB, const UShort& W, const double& deform )
{
	constexpr double	sqrt3by2 { 0.866025403784439 };
	double				dx, dy, stress {};

	for ( int i = 0; i < bonds.total; ++i )
	{
		dx = points.X[bonds.B[i]] + shiftB.bx[i] - points.X[bonds.A[i]] - shiftB.ax[i];
		dy = points.Y[bonds.B[i]] + shiftB.by[i] - points.Y[bonds.A[i]] - shiftB.ay[i];

		stress += ( 1.0 - (bonds.L0[i] / std::sqrt( dx*dx + dy*dy )) ) * dy * dx;
	}

	stress /= (W * W * sqrt3by2 * (1.0 + deform));
	return stress;
}

/* ENERGY MINIMISATION ROUTINE 1 - STEEPEST DESCENT --------------------------------------------------------------------- */
double SD_spring( POINT& points, const BOND& bonds, const SHIFT_BOND& shiftB )
{
	constexpr double	dt0 {0.1};
	constexpr uint		loopMax {static_cast<uint>(1e3)};
	double				energy, energyOLD;

	// Initialise system
	energyOLD = updateForce_spring( points, bonds, shiftB );

	if ( energyOLD < 1e-14 ) return energyOLD;

	for ( int loopNum = 0; loopNum < loopMax; ++loopNum )
	{
		// Explicit Euler Integration
		for ( int i = 0; i < points.total; ++i )
		{
			points.X[i] += points.FX[i] * dt0;
			points.Y[i] += points.FY[i] * dt0;
		}

		// Update current bond lengths & angles, calculate force on each node
		energy = updateForce_spring( points, bonds, shiftB );

		// Stopping condition
		//if ( std::abs(energyOLD - energy) < 1e-12 )	break;
		//else										energyOLD = energy;

		//if ( maxForce( points ) )	{std::cout << "MAXFORCE:\t" << i << '\n'; break;}
	}

	return energy;
}

/* RESCALE VELOCITIES USING ALPHA & INTEGRATE -------------------------------------------------------------------------------- */
void inline integrate_spring ( POINT& points, const double& alpha, const double& dt )
{
	const double	oneminusalpha = 1.0 - alpha;
	double			modVsq, modFsq, forceScale;

	for ( int i = 0; i < points.total; ++i )
	{
		// Add accelaration due to forces
		points.VX[i] += points.FX[i] * dt;
		points.VY[i] += points.FY[i] * dt;

		// Rescale - FIRE STEERING
 		modVsq = points.VX[i] * points.VX[i] + points.VY[i] * points.VY[i];
		modFsq = points.FX[i] * points.FX[i] + points.FY[i] * points.FY[i];

		if ( modFsq )
		{
			forceScale = alpha * std::sqrt( modVsq / modFsq );
			points.VX[i] = oneminusalpha * points.VX[i] + forceScale * points.FX[i];
			points.VY[i] = oneminusalpha * points.VY[i] + forceScale * points.FY[i];
		}

		// Integrate with rescaled velocity
		points.X[i] += points.VX[i] * dt;
		points.Y[i] += points.VY[i] * dt;
	}
}

/* RESET VELOCITIES ACCORDING TO FORCE -------------------------------------------------------------------------------------- */
void inline resetVelocity_spring ( POINT& points, const double& dt )
{
	for ( int i = 0; i < points.total; ++i )
	{
		points.VX[i] = points.FX[i] * dt;
		points.VY[i] = points.FY[i] * dt;
	}
}

/* ENERGY MINIMISATION ROUTINE 2 - FIRE ----------------------------------------------------------------------------------- */
double FIRE_spring(POINT& points, const BOND& bonds, const double& dt0, const UShort& W, SHIFT_BOND& shiftB, const double& strain, const double& deform)
{
	constexpr UShort	Ndelay {5};
	constexpr uint		loopMAX {static_cast<uint>(1e6)};
	constexpr double	fdec {0.5}, finc {1.1}, falp {0.99}, alpha0 {0.15};
	const double		dtmin {dt0*0.1}, dtmax{dt0*10.0};
	uint				loopNum {}, Npos {};
	double				energy, energyOLD, energyMIN, alpha {alpha0}, power, maxVel, dt {dt0};
	POINT				pointsMIN = points;

	// Random number generator
	//std::mt19937_64 						generator(std::random_device{}());
	//std::uniform_real_distribution<double>	distrib (-1.0, 1.0);

	// Set up mirror points
	shiftPointsB( bonds, shiftB, W, strain, deform );

	// Energy & Forces before minimising
	energyOLD = updateForce_spring( points, bonds, shiftB );
	energyMIN = energyOLD;
	if ( energyMIN < 1e-14 )	return energyMIN;	// No need to minimise

	// Initialise velocities according to force
	resetVelocity_spring( points, dt );

	// Standard FIRE step
	while ( loopNum < loopMAX )
	{
		power = 0.0;
		for ( int i=0; i<points.total; ++i ) power += points.FX[i] * points.VX[i] + points.FY[i] * points.VY[i];

		// If Power > 0 for more than Ndelay steps, increase dt & alpha
		if (power > 0.0)
		{
			Npos++;
			if ( Npos > Ndelay )
			{
				dt = std::min( dt*finc, dtmax );
				alpha *= falp;
			}
		}
		else	// Go half step back - reset force/velocity, reduce dt & alpha
		{
			dt = std::max( dt*fdec, dtmin );
			for ( int i = 0; i < points.total; ++i )
			{
				points.X[i] -= points.VX[i]*dt;
				points.Y[i]	-= points.VY[i]*dt;
			}

			Npos = 0;
			alpha = alpha0;
			energyOLD = updateForce_spring( points, bonds, shiftB );
			resetVelocity_spring( points, dt );
		}

		// Integration step
		integrate_spring( points, alpha, dt );

		energy = updateForce_spring( points, bonds, shiftB );

		// Continue along gradient
		if ( energy < energyOLD )
		{
			loopNum++;
			energyOLD = energy;

			if ( energy < energyMIN )
			{
				energyMIN = energy;
				pointsMIN = points;	// save minimum energy config
			}
		}
		else	// Overshoot - reduce dt, go half step back, & reset
 		{
 			if ( dt != dtmin )	dt = std::max(dt*fdec, dtmin);
			else				break;							// Exit

			for ( int i = 0; i < points.total; ++i )
			{
				points.X[i] -= points.VX[i] * dt;
				points.Y[i]	-= points.VY[i] * dt;
			}
	
			// Reset parameters
			energyOLD = updateForce_spring( points, bonds, shiftB );	// new baseline
			resetVelocity_spring( points, dt );
			loopNum++;
			Npos = 0;
			alpha = alpha0;
		}
	}
	// EXIT - return lowest energy reached at any point (not necessarily current) & set points to minima
	points = pointsMIN;
	return energyMIN;
}