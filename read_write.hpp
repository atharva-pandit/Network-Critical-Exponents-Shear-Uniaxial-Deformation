/* MAKE FILE EMPTY --------------------------------------------- */
void clearfile( const std::string& filename )
{
	std::ofstream outdata;
	outdata.open(filename);
	outdata.close();
}

/* CONVERT NUMBER TO STRING WITH FIXED PRECISION -------------------------------------------- */
template <typename T1>
std::string mkStr(const T1 num, const int p)
{
    std::ostringstream out;
    out.precision(p);
    out << std::fixed << num;
    return std::move(out).str();
}

/* GET NUMBER OF LINES IN A FILE ------------------------------------------- */
int filesize ( const std::string& filename )
{
	int	count {};
	std::string line;
    std::ifstream file( filename );

	if (file)	while ( std::getline(file, line) )	count++;
	else
	{
		std::cerr << "ERROR : " << filename << " not found !" << '\n';
		std::exit(1);
	}
	file.close();

	return count;
}

/* READ POINTS FROM FILE -------------------------------------------- */
void readPoint( const std::string& filename, POINT& points )
{
	std::ifstream file(filename);

	uint	id, count {};
	double	x,y;

	points.total = filesize( filename );
	points.allocate_memory();

	if (file) while (file >> id >> x >> y)
	{
		points.ID[count] = id;
		points.X[count] = x;		points.Y[count] = y;
		points.VX[count] = 0.0;		points.VY[count] = 0.0;
		points.FX[count] = 0.0;		points.FY[count] = 0.0;

		count++;
	}
	else
	{
		std::cerr << "ERROR : " << filename << " not found !" << '\n';
		std::exit(1);
	}
	file.close();
}

/* READ BONDS FROM FILE -------------------------------------------- */
void readBond( const std::string& filename, BOND& bonds )
{
	std::ifstream file(filename);

	UShort	prd;
	uint	id, a, b, count {};
	double	l0;

	bonds.total = filesize( filename );
	bonds.allocate_memory();

	if (file) while (file >> id >> a >> b >> l0 >> prd)
	{
			bonds.ID[count] = id;
			bonds.A[count] = a;
			bonds.B[count] = b;
			bonds.L0[count] = l0;
			bonds.PRD[count] = prd;
			count++;
	}
	else
	{
		std::cerr << "ERROR : " << filename << " not found !" << '\n';
		std::exit(1);
	}
	file.close();
}

/* READ ANGLES FROM FILE -------------------------------------------- */
void readAngle( const std::string& filename, ANGLE& angles )
{
	std::ifstream file( filename );

	UShort	prd1, prd2;
	uint	id, a, b, c, count {};
	double	t0;

	angles.total = filesize( filename );
	angles.allocate_memory();

	if (file) while (file >> id >> a >> b >> c >> t0 >> prd1 >> prd2)
	{
		angles.ID[count] = id;
		angles.A[count] = a;
		angles.B[count] = b;
		angles.C[count] = c;
		angles.T0[count] = t0;
		angles.PRD1[count] = prd1;
		angles.PRD2[count] = prd2;
		count++;
	}
	else
	{
		std::cerr << "ERROR : " << filename << " not found !" << '\n';
		std::exit(1);
	}
	file.close();
}

/* LAMMPS INPUT FILE -------------------------------------- (ADD=0 for ATOM / ADD=1 for POINT) */
template <typename T2, typename T3, typename T4>
void writeLammpsInp(const T2& points, const bool add, const T3& bonds, const T4& angles, double strain, double deform, UShort W, float& Z, UShort K, UShort& V, int t)
{
	// File Name -------------------------------------------------------------------
	const std::string lmpfile { "movies/W"+mkStr(W,0)+"_Z"+mkStr(Z,2)+"_K"+mkStr(K,0)+"_E"+mkStr(deform,2)+"_V"+mkStr(V,0)+"_t"+mkStr(t,0)+".dat" };
	std::ofstream od;
	od.open(lmpfile); // makes file empty, append with std::ios::app;
	constexpr double	sqrt3by2 { 0.866025403784439 };

	// Header ----------------------------------------------------------------
	od << "Triangular lattice under strain" << "\n\n";	// comment
	od << points.total << " atoms" << '\n';
	od << "1 atom types" << '\n';
	od << bonds.total << " bonds" << '\n';
	od << bonds.total << " bond types" << '\n';
	od << angles.total << " angles" << '\n';
	od << angles.total << " angle types" << "\n\n";	

	// Simulation box borders --------------------------------------------------
	od << -W*0.5*(1.0+sqrt3by2*(1.0+deform)*strain) << " " << W*0.5*(1.0-sqrt3by2*(1.0+deform)*strain) << " xlo xhi" << '\n';
	od << -W*0.5*(1.0+deform)*sqrt3by2 << " " << W*0.5*(1.0+deform)*sqrt3by2 << " ylo yhi" << '\n';
	od << "-0.5 0.5 zlo zhi" << '\n';
	od << W*sqrt3by2*(1.0L+deform)*strain << " 0.0 0.0 xy xz yz" << "\n\n";

	// Atomic Masses -----------------------------------------------------------------------
	od << "Masses\n\n";
	od << "1 1.00784" << "\n\n";

	// Write Bond & Angle Coeffs & rest values -----------------------------------------------
 	od << "Bond Coeffs # harmonic\n\n";									// bond-type  mu  rest-length
	for (int i=0; i<bonds.total; i++)	od << bonds.ID[i] << "\t1.0\t" << bonds.L0[i] << '\n';		// mu=1

	od << "\nAngle Coeffs # harmonic\n\n";								// angle-type  kappa  rest-angle
	for (int i=0; i<angles.total; i++)	od << angles.ID[i] << "\t0.001\t" << angles.T0[i] << '\n';	// kappa=0.001

	// Write Atom Positions ---------------------------------------------------
	od << "\nAtoms # angle\n\n";										// atom-id  molecule-id  atom-type  x  y  z
	for (int i=0; i<points.total; i++) od << points.ID[i] << '\t' << points.ID[i] << "\t1\t" << points.X[i] << '\t' << points.Y[i] << "\t0.0\n";

	// Write Bonds -------------------------------------------------------------
	od << "\nBonds\n\n";												// bond-id  bond-type  a  b
	for (int i=0; i<bonds.total; i++)	od << bonds.ID[i] << '\t' << bonds.ID[i] << '\t' << bonds.A[i]+add << '\t' << bonds.B[i]+add << '\n';

	// Write Angles ---------------------------------------------------------------
	od << "\nAngles\n\n";												// angle-id  angle-type  a  b  c
	for (int i=0; i<angles.total; i++) od << angles.ID[i] << '\t' << angles.ID[i] << '\t' << angles.A[i]+add << '\t' << angles.B[i]+add << '\t' << angles.C[i]+add << '\n';

	od.close();
}

/* LAMMPS STYLE OUTPUT DUMP FILE ------------------------------------------------------------------ */
void writeDump( POINT& p, int& tstep, double strain, double deform, const int& W, float& avgZ  )
{
	const std::string	filename {"p"+mkStr(W,0)+"_Z"+mkStr(avgZ,1)+"_t"+std::to_string(tstep)+".dat"};
	const double		sqrt3by2 {0.866025403784439}, xy { W*(0.5 + strain*sqrt3by2)*(1.0+deform) };

	std::ofstream	od;
	od.open(filename);

	od << "ITEM: TIMESTEP" << '\n';
	od << tstep << '\n';

	od << "ITEM: NUMBER OF ATOMS" << '\n';
	od << p.total << '\n';

	od << "ITEM: BOX BOUNDS xy xz yz pp pp pp" << '\n';
	od << "0.0 " << W + xy << " " << xy << '\n';				// xlo_bound xhi_bound xy
	od << "0.0 " << W*sqrt3by2*(1.0+deform) << " 0.0" << '\n';	// ylo_bound yhi_bound yz
	od << "-0.5 0.5 0.0" << '\n';								// zlo_bound zhi_bound yz

	od << "ITEM: ATOMS id x y vx vy fx fy" << '\n';
	for ( int i = 0; i < p.total; i++ )
	{
		od << p.ID[i] << " " << p.X[i] << " " << p.Y[i] << " " << p.VX[i] << " " << p.VY[i] << " " << p.FX[i] << " " << p.FY[i] << '\n';
	}
	
	od.close();
}

/* NORMAL INPUT FILE -------------------------------------------------------- */
void writeNormalInp( const ATOM& atoms, const BOND_NET& bonds, const ANGLE_NET& angles, UShort W, float& avgZ, UShort& V )
{
	const std::string pointfile	{"data/p"+mkStr(W,0)+"_Z"+mkStr(avgZ,2)+"_V"+mkStr(V,0)+".dat"},
					  bondfile	{"data/b"+mkStr(W,0)+"_Z"+mkStr(avgZ,2)+"_V"+mkStr(V,0)+".dat"},
					  anglefile	{"data/a"+mkStr(W,0)+"_Z"+mkStr(avgZ,2)+"_V"+mkStr(V,0)+".dat"};

	std::ofstream pf, bf, af;

	// Write positions to points file
	pf.open( pointfile );
	for (int i=0; i<atoms.total; i++) pf << std::fixed << std::setprecision(14) << atoms.ID[i] << '\t' << atoms.X[i] << '\t' << atoms.Y[i] << '\n';
	pf.close();

	// Write bonds to bonds file
	bf.open( bondfile );
	for (int i=0; i<bonds.total; i++) bf << std::fixed << std::setprecision(14) << bonds.ID[i] << '\t' << bonds.A[i]-1 << '\t' << bonds.B[i]-1 << '\t' << bonds.L0[i] << '\t' << bonds.PRD[i] << '\n';
	bf.close();

	// Write angles to angles file
	af.open( anglefile );
	for (int i=0; i<angles.total; i++)
	af << std::fixed << std::setprecision(14) << angles.ID[i] << '\t' << angles.A[i]-1 << '\t' << angles.B[i]-1 << '\t' << angles.C[i]-1 << '\t' << angles.T0[i] << '\t' << angles.PRD1[i] << '\t' << angles.PRD2[i] << '\n';
	af.close();
}

/* WRITE TIMESTEP & ENERGY ------------------------------------------------------------------ */
template <typename T3>
void writeEnergy( T3 step, double& value, const std::string& outfilename )
{
	std::ofstream od;
	od.open( outfilename, std::ios::app );

	od << std::fixed << std::setprecision(14) << step << '\t' << value << '\n';

	od.close();
}

/* WRITE TWO VALUES TO FILE ------------------------------------------------------------------ */
template <typename T4>
void writeTwo( T4 step, double& value1, double& value2, const std::string& outfilename )
{
	std::ofstream od;
	od.open( outfilename, std::ios::app );

	od << std::fixed << std::setprecision(14) << step << '\t' << value1 << '\t' << value2 << '\n';

	od.close();
}

/* SAVE CURRENT STATE ------------------------------------------------------------------ */
void saveState( const POINT& p, const UShort& W, float& Z, double& deform, UShort& K, UShort& V, double dirxn, int tstep )
{
	std::string savefile {"save/p"+mkStr(W,0)+"_Z"+mkStr(Z,2)+"_K"+mkStr(K,0)+"_E"+mkStr(deform,2)+"_V"+mkStr(V,0)+"_r"+mkStr(dirxn,0)+"_t"+mkStr(tstep,0)+".dat"};

	std::ofstream pf;
	pf.open( savefile );

	for (int i=0; i < p.total; i++)	pf << std::fixed << std::setprecision(14) << p.ID[i] << '\t' << p.X[i] << '\t' << p.Y[i] << '\n';

	pf.close();
}

void saveStateDef( const POINT& p, const UShort& W, float& Z, UShort& K, UShort& V, const double& deform )
{
	std::string savefile {"saveDEF/p"+mkStr(W,0)+"_Z"+mkStr(Z,2)+"_K"+mkStr(K,0)+"_V"+mkStr(V,0)+"_E"+mkStr(deform,3)+".dat"};

	std::ofstream pf;
	pf.open( savefile );

	for (int i=0; i < p.total; i++)	pf << std::fixed << std::setprecision(14) << p.ID[i] << '\t' << p.X[i] << '\t' << p.Y[i] << '\n';

	pf.close();
}