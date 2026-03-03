#include "setup.hpp"
#include "read_write.hpp"

int main ( int argc, char* argv[] )
{
	/* VARIABLES ---------------------------------------------------------------- */
	ATOM							atoms;
	BOND_NET						bonds;
	ANGLE_NET						angles;
	std::map<int,int>				network_count;
	std::map<int,int>::iterator		network_max;
	std::unordered_map<uint, uint>	idswap;
	std::array<uint,6>				bond;
	std::array<bool,4>				brtl, brtlA, brtlB;
	SHIFT_ANGLE						shiftA;
	SHIFT_BOND						shiftB;
	std::stringstream				convert1 {argv[1]}, convert2 {argv[2]}, convert3 {argv[3]};

	const double	sqrt3by2 { std::sqrt(3.0)*0.5 };
	double			x;
	bool			check, checkX, checkY;
	UShort			W, V, rowshift;
	uint			aid {}, bid {}, netID {}, ind, fil, temp1, temp2, idx;
	float			Z, avgZ {};
	Timer			Time;

	/* INPUT FROM TERMINAL ----------------------------------------------------------- */
	convert1 >> W;	convert2 >> Z;	convert3 >> V;
	int				total_atoms {W*W}, newID {total_atoms};
	const double	boxXmin {-W*0.5}, boxYmin {-(W-0.5)*0.5*sqrt3by2};

	/* RESERVE MEMORY ---------------------------------------------------------------*/
	reserve_memory( atoms, bonds, angles, total_atoms );

	/* CALCULATE LATTICE POSITIONS & BONDS --------------------------------------- */
	for (int j = 0; j < W; j++)
	{
		rowshift = std::abs( j % 2 );
		for (int i = 0; i < W; i++)
		{
			atoms.ID.push_back( i + j * W + 1 );
			atoms.X.push_back( boxXmin + i + 0.125 + 0.5*rowshift );
			atoms.Y.push_back( boxYmin + j * sqrt3by2 );
			atoms.Z.push_back( 0 );

			bond[0] = (i + rowshift)%W + ((j+1)*W)%total_atoms + 1;		// (i,j+1)
			bond[1] = (i+1)%W + j*W + 1;								// (i+1,j)
			bond[2] = (i + rowshift)%W + ((W+j-1)*W)%total_atoms + 1;	// (i,j-1)
			bond[3] = (W+i-1+rowshift)%W + ((W+j-1)*W)%total_atoms + 1;	// (i-1,j-1)
			bond[4] = (W+i-1)%W + j*W + 1;								// (i-1,j)
			bond[5] = (W+i-1+rowshift)%W + ((j+1)*W)%total_atoms + 1;	// (i-1,j+1)

			if ( j == 0 )		brtl[0] = 1;	else brtl[0] = 0;
			if ( i == (W-1) )	brtl[1] = 1;	else brtl[1] = 0;
			if ( j == (W-1) )	brtl[2] = 1;	else brtl[2] = 0;
			if ( i == 0 )		brtl[3] = 1;	else brtl[3] = 0;
	
			atoms.BOND.push_back( bond );
			atoms.BRTL.push_back( brtl );
			atoms.PHANTOM.push_back( 0 );
			atoms.NET.push_back(0);
		}
	}
	atoms.total = atoms.ID.size();

	/* ADD PHANTOM ATOMS AT EACH SITE ------------------------------------------------------------------ */
	std::mt19937_64 generator( std::random_device{}() );				// seed
	std::uniform_int_distribution<int> filament_distribution (0, 2);	// random int = 0,1,2
 
	for (int j=0; j<6; j++) bond[j] = 0;	// empty bond to be filled in

	for ( int i=0; i<total_atoms; i++ )
	{
		atoms.ID.push_back( ++newID );
		atoms.X.push_back( atoms.X[i] );
		atoms.Y.push_back( atoms.Y[i] );
		atoms.Z.push_back( 0 );

		for (int j=0; j<4; j++)	brtl[j] = atoms.BRTL[i][j];
		atoms.BOND.push_back( bond );	// empty
		atoms.BRTL.push_back( brtl );	// copied
		atoms.PHANTOM.push_back( 1 );
		atoms.NET.push_back(0);

		fil = filament_distribution( generator );	// select random filament to phantomise

		atoms.BOND[newID-1][fil] = atoms.BOND[i][fil];	// connect phantom to neighbours
		atoms.BOND[newID-1][fil+3] = atoms.BOND[i][fil+3];

		atoms.BOND[atoms.BOND[i][fil]-1][fil+3] = newID;	// connect neibs to phantom
		atoms.BOND[atoms.BOND[i][fil+3]-1][fil] = newID;

		atoms.BOND[i][fil] = 0;					// disconnect atom from neighbours
		atoms.BOND[i][fil+3] = 0;
	}
	atoms.total = atoms.ID.size();

	/* DELETE ONE BOND PER ROW & COLUMN - MAX FILAMENT LENGTH = SYSTEM SIZE ------------- */
 	std::uniform_int_distribution<int> distribution_row (0, W-1);

 	for (int i = 0; i < W; i++)
	{
		check = 1;
		while( check )
		{
			ind = i + W * distribution_row(generator);
			if ( atoms.BOND[ind][0] && !atoms.PHANTOM[atoms.BOND[ind][0]-1] )
			{
				atoms.BOND[atoms.BOND[ind][0]-1][3] = 0;
				atoms.BOND[ind][0] = 0;
				check = 0;
			}
		}
	}
	for (int j = 0; j < W; j++)
	{
		check = 1;
		while( check )
		{
			ind = j * W + distribution_row(generator);
			if ( atoms.BOND[ind][1] && !atoms.PHANTOM[atoms.BOND[ind][1]-1] )
			{
				atoms.BOND[atoms.BOND[ind][1]-1][4] = 0;
				atoms.BOND[ind][1] = 0;
				check = 0;
			}
		}
	}

	std::vector<int> reverse_diagonal_ind {0};
	for (int k = 1; k < W; k++) reverse_diagonal_ind.emplace_back(W + k*(W-1));
	for (int l = 0; l < W; l++)
	{
		check = 1;
		while( check )
		{
			ind = (l*W + reverse_diagonal_ind[distribution_row(generator)]) % total_atoms;
			if ( atoms.BOND[ind][2] && !atoms.PHANTOM[atoms.BOND[ind][2]-1] )
			{
				atoms.BOND[atoms.BOND[ind][2]-1][5] = 0;
				atoms.BOND[ind][2] = 0;
				check = 0;
			}
		}
	}

	/* DELETE RANDOM BONDS UNTIL DESIRED CONNECTIVITY ---------------------------- */
 	std::uniform_int_distribution<int> atom_distribution ( 0, atoms.total-1 );

	localZ( atoms );		// set Z for every atom
	avgZ = calcZ( atoms );	// network connectivity

	while ( avgZ > Z )
	{
		ind = atom_distribution( generator );
		fil = filament_distribution( generator );
		temp1 = atoms.BOND[ind][fil];

		if ( temp1 )
		{
			atoms.BOND[temp1-1][fil+3] = 0;	// disconnect from the other side first
			atoms.Z[temp1-1] -= 1;
			atoms.BOND[ind][fil] = 0;		// delete bond
			atoms.Z[ind] -= 1;

			// remove atoms if they are down to a single bond
			if ( atoms.Z[ind] == 1 )		remove_dangling_ends( ind, atoms );
			if ( atoms.Z[temp1-1] == 1 )	remove_dangling_ends( temp1-1, atoms );

			avgZ = calcZ( atoms );				// update network connectivity
		}
	}

	/* LABEL DISJOINT NETWORK SETS ------------------------------------------------ */
 	for ( int i = 0; i < atoms.total; i++ )
	{
		if ( atoms.NET[i] == 0 )
		{
			atoms.NET[i] = ++netID;
			network_count.insert( {netID, 0} );
			set_network ( i, atoms );
		}
		network_count[ atoms.NET[i] ]++;
	}
	network_max = std::max_element ( network_count.begin(), network_count.end(),
					[]( const std::map<int,int>::value_type& a, const std::map<int,int>::value_type& b )
					{ return a.second < b.second; } );
					
	/* DELETE NON-LARGEST SUB-NETWORKS ----------------------------------------------- */
	for ( int i = 0; i < atoms.total; i++ )
	{
		if ( atoms.NET[i] != network_max->first )	atoms.ID[i] = 0;
	}

	/* DELETE VOID ATOMS FROM LIST ------------------------------------------------------ */
	idx = 0;	std::erase_if( atoms.X, [&](double) { return atoms.ID[idx++] == 0; } );
	idx = 0;	std::erase_if( atoms.Y, [&](double) { return atoms.ID[idx++] == 0; } );
	idx = 0;	std::erase_if( atoms.Z, [&](float)	{ return atoms.ID[idx++] == 0; } );
	idx = 0;	std::erase_if( atoms.PHANTOM, [&](bool) { return atoms.ID[idx++] == 0; } );
	idx = 0;	std::erase_if( atoms.BOND, [&](std::array<uint,6>) { return atoms.ID[idx++] == 0; } );
	idx = 0;	std::erase_if( atoms.BRTL, [&](std::array<bool,4>) { return atoms.ID[idx++] == 0; } );
	idx = 0;	std::erase_if( atoms.NET, [&](uint) { return atoms.ID[idx++] == 0; } );
	idx = 0;	std::erase_if( atoms.ID, [&](uint)	{ return atoms.ID[idx++] == 0; } );
	atoms.total = atoms.ID.size();

	/* REASSIGN NEW INDICES TO ATOMS & BONDS ----------------------------------------------------- */
	idswap.insert( {0, 0} );
 	for ( int i = 0; i < atoms.total; i++)
	{
		idswap.insert( {atoms.ID[i], i+1} );
		atoms.ID[i] = i+1;
	}

	for ( int i = 0; i < atoms.total; i++)
	{
		for ( int j = 0; j < 6; j++ )	if ( atoms.BOND[i][j] ) atoms.BOND[i][j] = idswap[ atoms.BOND[i][j] ];
	}

	// Final connectivity
	localZ( atoms );
	avgZ = calcZ( atoms );

	/* ATOM IN CENTRE OF BONDS FOR BUCKLING ----------------------------------------------------- */
   	newID = atoms.total;
	for (int k=0; k<4; k++)	brtl[k] = 0; 
	for (int k=0; k<6; k++)	bond[k] = 0; 

	for ( int i = 0; i < atoms.total; i++ )
	{
		brtlA = atoms.BRTL[i];

		for ( int j = 0; j < 3; j++ )
		{
			if ( atoms.BOND[i][j] )
			{
				temp1 = atoms.BOND[i][j] - 1;
				brtlB = atoms.BRTL[ temp1 ];

				atoms.ID.push_back( ++newID );
				atoms.PHANTOM.push_back ( 0 );
				atoms.BRTL.push_back( brtl );	// empty
				checkX = 1; checkY = 1;

				if ( (brtlA[1] && brtlB[3]) || (brtlA[3] && brtlB[1]) )
				{
					x = 0.5 * ( atoms.X[i] + atoms.X[temp1] + W );
					if ( x > W*0.5 ) x -= W;
					atoms.X.push_back( x );
					checkX = 0;
				}
				if ( (brtlA[2] && brtlB[0]) || (brtlA[0] && brtlB[2]) )
				{
					atoms.Y.push_back( 0.5*( atoms.Y[i] + atoms.Y[temp1] + W*sqrt3by2 ) );
					checkY = 0;
				}
				if (checkX)	atoms.X.push_back( 0.5 * ( atoms.X[i] + atoms.X[temp1] ) );
				if (checkY)	atoms.Y.push_back( 0.5 * ( atoms.Y[i] + atoms.Y[temp1] ) );

				atoms.BOND.push_back( bond );
				atoms.BOND[temp1][j+3] = newID;
				atoms.BOND[i][j] = newID;
				atoms.BOND[newID-1][j+3] = atoms.ID[i];
				atoms.BOND[newID-1][j] = atoms.ID[temp1];
			}
		}
	}
	atoms.total = atoms.ID.size();

	/* RESET ALL ATOM BOUNDARY LABELS ------------------------------------------------- */
	for ( std::array<bool,4>& arr : atoms.BRTL ) for ( bool& val : arr ) val = 0; // set all values to 0

	for ( int i = 0; i < atoms.total; i++ )
	{
		if ( atoms.X[i] <= -W*0.5 + 0.25 )					 atoms.BRTL[i][3] = 1;
		else if ( atoms.X[i] >= W*0.5 - 0.5 )				 atoms.BRTL[i][1] = 1;
		if ( atoms.Y[i] < -(W-0.5)*0.5*sqrt3by2 + 0.25 )	 atoms.BRTL[i][0] = 1;
		else if ( atoms.Y[i] > (W-0.5)*0.5*sqrt3by2 - 0.25 ) atoms.BRTL[i][2] = 1;
	}

	/* CREATE BONDS LIST ---------------------------------------------- */
	for ( int i = 0; i < atoms.total; i++ )
	{
		for ( int j = 0; j < 3; j++)
		{
			temp1 = atoms.BOND[i][j];
			if ( temp1 )
			{
				bonds.ID.push_back( ++bid );
				bonds.A.push_back( atoms.ID[i] );
				bonds.B.push_back( temp1 );
				bonds.PRD.push_back( periodic( atoms.X[i], atoms.Y[i], atoms.X[temp1-1], atoms.Y[temp1-1], atoms.BRTL[i], atoms.BRTL[temp1-1], W ) );
			}
		}
	}
	bonds.total = bonds.ID.size();

	/* CREATE ANGLE LIST ------------------------------------------------------- */
	for ( int i = 0; i < atoms.total; i++ )
	{
		for ( int j = 0; j < 3; j++)
		{
			temp1 = atoms.BOND[i][j];
			temp2 = atoms.BOND[i][j+3];
			if ( temp1 && temp2 )
			{
				angles.ID.push_back( ++aid );
				angles.A.push_back( temp1 );
				angles.B.push_back( atoms.ID[i] );
				angles.C.push_back( temp2 );
				angles.PRD1.push_back( periodic( atoms.X[temp1-1], atoms.Y[temp1-1], atoms.X[i], atoms.Y[i], atoms.BRTL[temp1-1], atoms.BRTL[i], W ) );
				angles.PRD2.push_back( periodic( atoms.X[i], atoms.Y[i], atoms.X[temp2-1], atoms.Y[temp2-1], atoms.BRTL[i], atoms.BRTL[temp2-1], W ) );
			}
		}
	}
	angles.total = angles.ID.size();

	/* CALCULATE REST VALUES FOR BONDS & ANGLES --------------------------------------- */
	shiftA.allocate_memory( angles.total );
	shiftB.allocate_memory( bonds.total );
	shiftPointsA( angles, shiftA, W, 0.0, 0.0 );
	shiftPointsB( bonds, shiftB, W, 0.0, 0.0 );
	nDist( atoms, bonds, shiftB, W );		// zero strain & extension
	nDegrees( atoms, angles, shiftA, W );

	/* SAVE NETWORK ----------------------------------------------------------------- */
	writeNormalInp( atoms, bonds, angles, W, avgZ, V );
	//writeLammpsInp( atoms, 0, bonds, angles, 0.0, 0.0, W, avgZ, 0, V, 0 );

	/* DONE ---------------------------------------------------------------------- */
	std::cout << "W: " << W << "\tZ: "  << avgZ << "\tV: " << V << '\n';
	std::cout << "#Atoms : " << atoms.total << '\t' << "(#Deleted : " << 5*W*W - atoms.total << ")\n";
	std::cout << "#Bonds : " << bid << '\t' << "(#Deleted : " << 6*W*W - bid << ")\n";
	std::cout << "#Angles : " << aid << '\t' << "(#Deleted : " << 6*W*W - aid << ")\n";
	std::cout << "Time elapsed: " << Time.elapsed() << " seconds.\n";
	return 0;
}