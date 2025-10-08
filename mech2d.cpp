#include "inmost.h"

#define DAT double 

using namespace INMOST;

// Assuming a uniform Cartesian grid
// 
//          node(i,j+1)      face_y(i,j+1)  node(i+1,j+1)
//                      *--------+--------*
//                      |                 |
//                      |                 |
//                      |                 |
//          face_x(i,j) +     cell(i,j)   + face_x(i+1,j)
//                      |                 |
//                      |                 |
//                      |                 |
//                      *--------+--------*
//              node(i,j)    face_y(i,j)    node(i+1,j)
//
//
//   Therefore, for a NxN square cells grid,
//   Cell   index range: 0...N-1, 0...N-1, 0...N-1
//   Node   index range: 0...N,   0...N,   0...N
//   Face-x index range: 0...N,   0...N-1, 0...N-1
//   Face-y index range: 0...N-1, 0...N,   0...N-1
//   Face-z index range: 0...N-1, 0...N-1, 0...N
//
//   Cell   vars: sigma_xx, sigma_yy, sigma_zz
//   Node   vars: sigma_xy, sigma_xz, sigma_yz
//   Face-x vars: u_x
//   Face-y vars: u_y
//   Face-z vars: u_z

class Problem
{
private:
	// =====================================================================================================
	// Problem parameters: discretization and physical parameters
	int N;
	double dx;
	double dy;
	double dz;
	const double Lx = 1.0;
	const double Ly = Lx;
	const double Lz = Lx;
	const double E = 3.5e6;                               // Young's modulus, MPa
	const double nu = 0.3;                                // Poisson ratio
	const double lam = E * nu / (1 + nu) / (1 - 2 * nu);  // Lame parameter lambda
	const double mu = E / 2 / (1 + nu);                   // Lame parameter mu


	// =====================================================================================================
	// INMOST functionality
	/// Solution vector that contains all the unknowns: 
	/// sigma_xx, sigma_yy, sigma_zz, sigma_xy, sigma_xz, sigma_yz, sigma_yx, sigma_zx, sigma_zy, u_x, u_y, u_z
	INMOST::Sparse::Vector sol;
	/// Linear solver for systems with Jacobian matrix
	//INMOST::Solver S;
	/// Nonlinear residual
	//INMOST::Residual R;


	// =====================================================================================================
	// Index functions that map (i,j) for a specific unknown to the position in the global vector 'sol'
	int Isxx(int i, int j, int k); // sigma_xx
	int Isyy(int i, int j, int k); // sigma_yy
	int Iszz(int i, int j, int k); // sigma_zz
	int Isxy(int i, int j, int k); // sigma_xy
	int Isxz(int i, int j, int k); // sigma_xz
	int Isyz(int i, int j, int k); // sigma_yz
	int Isyx(int i, int j, int k); // sigma_yx
	int Iszx(int i, int j, int k); // sigma_zx
	int Iszy(int i, int j, int k); // sigma_zy
	int Iux(int i, int j, int k);  // u_x
	int Iuy(int i, int j, int k);  // u_y
	int Iuz(int i, int j, int k);  // u_z


	// =====================================================================================================
	// Functions that construct unknonwns from given locations in 'sol'
	unknown sxx(int i, int j, int k) { return unknown(sol[Isxx(i,j,k)], Isxx(i,j,k)); }
	unknown syy(int i, int j, int k) { return unknown(sol[Isyy(i,j,k)], Isyy(i,j,k)); }
	unknown szz(int i, int j, int k) { return unknown(sol[Iszz(i,j,k)], Iszz(i,j,k)); }
	unknown sxy(int i, int j, int k) { return unknown(sol[Isxy(i,j,k)], Isxy(i,j,k)); }
	unknown sxz(int i, int j, int k) { return unknown(sol[Isxz(i,j,k)], Isxz(i,j,k)); }
	unknown syz(int i, int j, int k) { return unknown(sol[Isyz(i,j,k)], Isyz(i,j,k)); }
	unknown syx(int i, int j, int k) { return unknown(sol[Isyx(i,j,k)], Isyx(i,j,k)); }
	unknown szx(int i, int j, int k) { return unknown(sol[Iszx(i,j,k)], Iszx(i,j,k)); }
	unknown szy(int i, int j, int k) { return unknown(sol[Iszy(i,j,k)], Iszy(i,j,k)); }
	unknown ux(int i, int j, int k)  { return unknown(sol[Iux(i,j,k)],  Iux(i,j,k));  }
	unknown uy(int i, int j, int k)  { return unknown(sol[Iuy(i,j,k)],  Iuy(i,j,k));  }
	unknown uz(int i, int j, int k)  { return unknown(sol[Iuz(i,j,k)],  Iuz(i,j,k));  }



public:
	Problem(int N_)
	{
		N = N_;
		dx = Lx / N;
		dy = Ly / N;
		dz = Lz / N;
	}
	~Problem() {}
	void run();
	void fillResidual(Residual &R);
	void saveVTK2D();
	void saveVTK3D();
};


int Problem::Isxx(int i, int j, int k)
{
	if (i < 0 || i > N-1) 
		std::cout << "Isxx: wrong i = " << i << " for N = " << N << std::endl;
	if (j < 0 || j > N-1)
		std::cout << "Isxx: wrong j = " << j << " for N = " << N << std::endl;
	if (k < 0 || k > N-1)
		std::cout << "Isxx: wrong k = " << k << " for N = " << N << std::endl;
	return (i*N + j)*N + k;
}

int Problem::Isyy(int i, int j, int k)
{
	if (i < 0 || i > N-1)
		std::cout << "Isyy: wrong i = " << i << " for N = " << N << std::endl;
	if (j < 0 || j > N-1)
		std::cout << "Isyy: wrong j = " << j << " for N = " << N << std::endl;
	if (k < 0 || k > N-1)
		std::cout << "Isyy: wrong k = " << k << " for N = " << N << std::endl;
	return (i*N + j)*N + k + N*N*N;
}

int Problem::Iszz(int i, int j, int k)
{
	if (i < 0 || i > N-1)
		std::cout << "Iszz: wrong i = " << i << " for N = " << N << std::endl;
	if (j < 0 || j > N-1)
		std::cout << "Iszz: wrong j = " << j << " for N = " << N << std::endl;
	if (k < 0 || k > N-1)
		std::cout << "Iszz: wrong k = " << k << " for N = " << N << std::endl;
	return (i*N + j)*N + k + N*N*N*2;
}

int Problem::Isxy(int i, int j, int k)
{
	if (i < 0 || i > N)
		std::cout << "Isxy: wrong i = " << i << " for N = " << N << std::endl;
	if (j < 0 || j > N)
		std::cout << "Isxy: wrong j = " << j << " for N = " << N << std::endl;
	if (k < 0 || k > N)
		std::cout << "Isxy: wrong k = " << k << " for N = " << N << std::endl;
	return (i*(N+1) + j)*(N+1) + k + N*N*N*3;
}

int Problem::Isxz(int i, int j, int k)
{
	if (i < 0 || i > N)
		std::cout << "Isxz: wrong i = " << i << " for N = " << N << std::endl;
	if (j < 0 || j > N)
		std::cout << "Isxz: wrong j = " << j << " for N = " << N << std::endl;
	if (k < 0 || k > N)
		std::cout << "Isxz: wrong k = " << k << " for N = " << N << std::endl;
	return (i*(N+1) + j)*(N+1) + k + N*N*N*3 + (N+1)*(N+1)*(N+1);
}

int Problem::Isyz(int i, int j, int k)
{
	if (i < 0 || i > N)
		std::cout << "Isyz: wrong i = " << i << " for N = " << N << std::endl;
	if (j < 0 || j > N)
		std::cout << "Isyz: wrong j = " << j << " for N = " << N << std::endl;
	if (k < 0 || k > N)
		std::cout << "Isyz: wrong k = " << k << " for N = " << N << std::endl;
	return (i*(N+1) + j)*(N+1) + k + N*N*N*3 + (N+1)*(N+1)*(N+1)*2;
}

int Problem::Isyx(int i, int j, int k)
{
	if (i < 0 || i > N)
		std::cout << "Isyx: wrong i = " << i << " for N = " << N << std::endl;
	if (j < 0 || j > N)
		std::cout << "Isyx: wrong j = " << j << " for N = " << N << std::endl;
	if (k < 0 || k > N)
		std::cout << "Isyx: wrong k = " << k << " for N = " << N << std::endl;
	return (i*(N+1) + j)*(N+1) + k + N*N*N*3 + (N+1)*(N+1)*(N+1)*3;
}

int Problem::Iszx(int i, int j, int k)
{
	if (i < 0 || i > N)
		std::cout << "Iszx: wrong i = " << i << " for N = " << N << std::endl;
	if (j < 0 || j > N)
		std::cout << "Iszx: wrong j = " << j << " for N = " << N << std::endl;
	if (k < 0 || k > N)
		std::cout << "Iszx: wrong k = " << k << " for N = " << N << std::endl;
	return (i*(N+1) + j)*(N+1) + k + N*N*N*3 + (N+1)*(N+1)*(N+1)*4;
}

int Problem::Iszy(int i, int j, int k)
{
	if (i < 0 || i > N)
		std::cout << "Iszy: wrong i = " << i << " for N = " << N << std::endl;
	if (j < 0 || j > N)
		std::cout << "Iszy: wrong j = " << j << " for N = " << N << std::endl;
	if (k < 0 || k > N)
		std::cout << "Iszy: wrong k = " << k << " for N = " << N << std::endl;
	return (i*(N+1) + j)*(N+1) + k + N*N*N*3 + (N+1)*(N+1)*(N+1)*5;
}

int Problem::Iux(int i, int j, int k)
{
	if (i < 0 || i > N)
		std::cout << "Iux: wrong i = " << i << " for N = " << N << std::endl;
	if (j < 0 || j > N-1)
		std::cout << "Iux: wrong j = " << j << " for N = " << N << std::endl;
	if (k < 0 || k > N-1)
		std::cout << "Iux: wrong k = " << k << " for N = " << N << std::endl;
	return (i*N + j)*(N+1) + k + N*N*N*3 + (N+1)*(N+1)*(N+1)*6;
}

int Problem::Iuy(int i, int j, int k)
{
	if (i < 0 || i > N-1)
		std::cout << "Iuy: wrong i = " << i << " for N = " << N << std::endl;
	if (j < 0 || j > N)
		std::cout << "Iuy: wrong j = " << j << " for N = " << N << std::endl;
	if (k < 0 || k > N-1)
		std::cout << "Iuy: wrong k = " << k << " for N = " << N << std::endl;
	return (i*(N+1) + j)*N + k + N*N*N*3 + (N+1)*(N+1)*(N+1)*6 + N*N*(N+1);
}

int Problem::Iuz(int i, int j, int k)
{
	if (i < 0 || i > N-1)
		std::cout << "Iuz: wrong i = " << i << " for N = " << N << std::endl;
	if (j < 0 || j > N-1)
		std::cout << "Iuz: wrong j = " << j << " for N = " << N << std::endl;
	if (k < 0 || k > N)
		std::cout << "Iuz: wrong k = " << k << " for N = " << N << std::endl;
	return (i*N + j)*N + k + N*N*N*3 + (N+1)*(N+1)*(N+1)*6 + N*N*(N+1)*2;
}

void Problem::saveVTK2D()
{
	
	// ===================================================================================
	// Save VTK with results (legacy code)
	std::ofstream out;
	out.open("res2d.vtk");
	if (!out.is_open()) {
		std::cout << "Couldn't open file res.vtk!";
		exit(0);
	}
	int nx = N, ny = N;
	out << "# vtk DataFile Version 3.0" << std::endl << std::endl;
	out << "ASCII" << std::endl;
	out << "DATASET STRUCTURED_GRID" << std::endl;
	out << "DIMENSIONS " << nx + 1 << " " << ny + 1 << " 1" << std::endl;
	out << "POINTS " << (nx + 1) * (ny + 1) << " DOUBLE" << std::endl;
	for (int j = 0; j <= ny; j++) {
		for (int i = 0; i <= nx; i++) {
			double uxij = 0.0, uxijm1 = 0.0;
			double uyij = 0.0, uyim1j = 0.0;
			if (i < nx && j < ny)
				uxij = sol[Iux(i, j, 3)];
			if (i < nx && j > 0)
				uxijm1 = sol[Iux(i, j - 1, 3)];
			if (i < nx && j < ny)
				uyij = sol[Iuy(i, j, 3)];
			if (i > 0 && j < ny)
				uyim1j = sol[Iuy(i - 1, j, 3)];

			double a = 0.0;
			double ux = 0.5 * (uxij + uxijm1);
			double uy = 0.5 * (uyij + uyim1j);
			out << i * dx + uy*a << " " << j * dy + uy*a<< " 0.0" << std::endl;
		}
	}

	out << "CELL_DATA " << nx * ny << std::endl;


	out << "SCALARS Stress_yy double" << std::endl;
	out << "LOOKUP_TABLE default" << std::endl;
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			out << sol[Isyy(i,j, 3)] << std::endl;
		}
	}

	out << "SCALARS Stress_xx double" << std::endl;
	out << "LOOKUP_TABLE default" << std::endl;
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			out << sol[Isxx(i,j, 3)] << std::endl;
		}
	}

	// Values in cell are averaged from cell nodes
	out << "SCALARS Txy double" << std::endl;
	out << "LOOKUP_TABLE default" << std::endl;
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			//out << 0.25 * (Txy[i + j * (nx + 1)] + Txy[i + j * (nx + 1)]
			//	+ Txy[i + 1 + j * (nx + 1)] + Txy[i + (j + 1) * (nx + 1)]) << std::endl;
			out << 0.25 * (sol[Isxy(i,j, 3)] + sol[Isxy(i+1, j, 3)]
				+ sol[Isxy(i, j+1, 3)] + sol[Isxy(i+1, j+1, 3)]) << std::endl;
		}
	}

	out << "POINT_DATA " << (nx+1) * (ny+1) << std::endl;

	out << "VECTORS Displacement double" << std::endl;
	for (int j = 0; j < ny+1; j++) {
		for (int i = 0; i < nx+1; i++) {
    
//          face_x(i,j) +     
//                      |       
//                      |       
//                      |         
//            -+--------*--------+-
//  face_y(i-1,j)   node(i,j)    face_y(i,j) 
//						|
//						|
//						+
//                    face_x(i,j-1)
			double uxij = 0.0, uxijm1 = 0.0;
			double uyij = 0.0, uyim1j = 0.0;
			if (i < nx && j < ny)
				uxij = sol[Iux(i,j, 3)]; 
			if (i < nx && j > 0)
				uxijm1 = sol[Iux(i,j-1, 3)];
			if (i < nx && j < ny)
				uyij = sol[Iuy(i,j, 3)];
			if (i > 0 && j < ny)
				uyim1j = sol[Iuy(i-1,j, 3)];

			double ux = 0.5 * (uxij + uxijm1);
			double uy = 0.5 * (uyij + uyim1j);

			out <<  ux << " " << uy << " 0.0" << std::endl;
		}
	}

	out.close();
}

void Problem::saveVTK3D()
{
	
	// ===================================================================================
	// Save VTK with results (legacy code)
	std::ofstream out;
	out.open("res3d.vtk");
	if (!out.is_open()) {
		std::cout << "Couldn't open file res.vtk!";
		exit(0);
	}
	int nx = N, ny = N, nz = N;
	out << "# vtk DataFile Version 3.0" << std::endl << std::endl;
	out << "ASCII" << std::endl;
	out << "DATASET STRUCTURED_GRID" << std::endl;
	out << "DIMENSIONS " << nx + 1 << " " << ny + 1 << " " << nz + 1 << std::endl;
	out << "POINTS " << (nx + 1) * (ny + 1) * (nz + 1) << " DOUBLE" << std::endl;
	for (int k = 0; k <= nz; k++)
		for (int j = 0; j <= ny; j++) 
			for (int i = 0; i <= nx; i++) 
				out << i * dx  << " " << j * dy << " " << k * dz << std::endl;

	out << "CELL_DATA " << nx * ny * nz << std::endl;


	out << "SCALARS Stress_xx double" << std::endl;
	out << "LOOKUP_TABLE default" << std::endl;
	for (int k = 0; k < nz; k++)
		for (int j = 0; j < ny; j++) 
			for (int i = 0; i < nx; i++) 
				out << sol[Isxx(i,j,k)] << std::endl;

	out << "SCALARS Stress_yy double" << std::endl;
	out << "LOOKUP_TABLE default" << std::endl;
	
	for (int k = 0; k < nz; k++)
		for (int j = 0; j < ny; j++)
			for (int i = 0; i < nx; i++)
				out << sol[Isyy(i,j,k)] << std::endl;

	out << "SCALARS Stress_zz double" << std::endl;
	out << "LOOKUP_TABLE default" << std::endl;
	for (int k = 0; k < nz; k++)
		for (int j = 0; j < ny; j++)
			for (int i = 0; i < nx; i++)
				out << sol[Iszz(i,j,k)] << std::endl;

	// Values in cell are averaged from cell nodes
	/*out << "SCALARS Stress_xy double" << std::endl;
	out << "LOOKUP_TABLE default" << std::endl;
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			//out << 0.25 * (Txy[i + j * (nx + 1)] + Txy[i + j * (nx + 1)]
			//	+ Txy[i + 1 + j * (nx + 1)] + Txy[i + (j + 1) * (nx + 1)]) << std::endl;
			out << 0.25 * (sol[Isxy(i,j, 3)] + sol[Isxy(i+1, j, 3)]
				+ sol[Isxy(i, j+1, 3)] + sol[Isxy(i+1, j+1, 3)]) << std::endl;
		}
	}*/

	out << "POINT_DATA " << (nx+1) * (ny+1) * (nz+1) << std::endl;

	out << "VECTORS Displacement double" << std::endl;
	for (int k = 0; k < nz+1; k++) {
		for (int j = 0; j < ny+1; j++) {
			for (int i = 0; i < nx+1; i++){
				//          face_x(i,j) +     
				//                      |       
				//                      |       
				//                      |         
				//            -+--------*--------+-
				//  face_y(i-1,j)   node(i,j)    face_y(i,j) 
				//						|
				//						|
				//						+
				//                    face_x(i,j-1)

				// u_x: averaging over 4 faces
				double uxijk = 0.0, uxijm1k = 0.0, uxijkm1 = 0.0, uxijm1km1 = 0.0;
				if (i < nx && j < ny && k < nz)
					uxijk = sol[Iux(i,j,k)]; 
				if (i < nx && j > 0 && k < nz)
					uxijm1k = sol[Iux(i,j-1,k)];
				if (i < nx && j < nx && k > 0)
					uxijkm1 = sol[Iux(i,j,k-1)];
				if (i < nx && j > 0 && k > 0)
					uxijm1km1 = sol[Iux(i,j-1,k-1)];

					
				double uyijk = 0.0, uyim1jk = 0.0, uyijkm1 = 0.0, uyim1jkm1 = 0.0;
				if (i < nx && j < ny && k < nz)
					uyijk = sol[Iuy(i,j,k)];
				if (i > 0 && j < ny && k < nz)
					uyim1jk = sol[Iuy(i-1,j,k)];
				if (i < nx && j < ny && k >0)
					uyijkm1 = sol[Iuy(i,j,k-1)];
				if (i > 0 && j < ny && k > 0)
					uyim1jkm1 = sol[Iuy(i-1,j,k-1)];

				double uzijk = 0.0, uzim1jk = 0.0, uzijm1k = 0.0, uzim1jm1k = 0.0;
				if (i < nx && j < ny && k < nz)
					uzijk = sol[Iuz(i,j,k)]; 
				if (i > 0 && j < ny && k < nz)
					uzim1jk = sol[Iuz(i-1,j,k)];
				if (i < nx && j > 0 && k < nz)
					uzijm1k = sol[Iuz(i,j-1,k)];
				if (i > 0 && j > 0 && k < nz)
					uzim1jm1k = sol[Iuz(i-1,j-1,k)];	

				double ux = 0.25 * (uxijk + uxijm1k + uxijkm1 + uxijm1km1);
				double uy = 0.25 * (uyijk + uyim1jk + uyijkm1 + uyim1jkm1);
				double uz = 0.25 * (uzijk + uzim1jk + uzijm1k + uzim1jm1k);

				out <<  ux << " " << uy << " " << uz << std::endl;
			}
		}
	}

	out << "SCALARS Stress_xy double" << std::endl;
	out << "LOOKUP_TABLE default" << std::endl;
	for (int k = 0; k <= nz; k++)
		for (int j = 0; j <= ny; j++)
			for (int i = 0; i <= nx; i++) 
				out << sol[Isxy(i, j, k)] << std::endl;

	out << "SCALARS Stress_xz double" << std::endl;
	out << "LOOKUP_TABLE default" << std::endl;
	for (int k = 0; k <= nz; k++)
		for (int j = 0; j <= ny; j++)
			for (int i = 0; i <= nx; i++) 
				out << sol[Isxz(i, j, k)] << std::endl;

	
	out << "SCALARS Stress_yz double" << std::endl;
	out << "LOOKUP_TABLE default" << std::endl;
	for (int k = 0; k <= nz; k++)
		for (int j = 0; j <= ny; j++)
			for (int i = 0; i <= nx; i++) 
				out << sol[Isyz(i, j, k)] << std::endl;

	
	out << "SCALARS Stress_yx double" << std::endl;
	out << "LOOKUP_TABLE default" << std::endl;
	for (int k = 0; k <= nz; k++)
		for (int j = 0; j <= ny; j++)
			for (int i = 0; i <= nx; i++) 
				out << sol[Isyx(i, j, k)] << std::endl;

	
	out << "SCALARS Stress_zx double" << std::endl;
	out << "LOOKUP_TABLE default" << std::endl;
	for (int k = 0; k <= nz; k++)
		for (int j = 0; j <= ny; j++)
			for (int i = 0; i <= nx; i++) 
				out << sol[Iszx(i, j, k)] << std::endl;

	
	out << "SCALARS Stress_zy double" << std::endl;
	out << "LOOKUP_TABLE default" << std::endl;
	for (int k = 0; k <= nz; k++)
		for (int j = 0; j <= ny; j++)
			for (int i = 0; i <= nx; i++) 
				out << sol[Iszy(i, j, k)] << std::endl;

	out.close();
}

void Problem::fillResidual(Residual &R)
{
	R.Clear();
	// --------------------------------- Cell loop
	// Equations for sxx, syy, szz:
	// 
	// sigma = 2*mu*eps + lam*tr(eps)*I (diagonal parts)
	//
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for(int k = 0; k < N; k++){
				variable uxijk = 0.0, uxip1jk = 0.0;
				variable uyijk = 0.0, uyijp1k = 0.0;
				variable uzijk = 0.0, uzijkp1 = 0.0;
				if (i > 0)
					uxijk = ux(i, j, k);
				if (i < N+1)
					uxip1jk = ux(i+1, j, k);
				if (j > 0)
					uyijk = uy(i, j, k);
				if (j < N+1)
					uyijp1k = uy(i, j+1, k);
				if (k > 0)
					uzijk = uz(i, j, k);
				if (k < N+1)
					uzijkp1 = uz(i, j, k+1);
				variable duxdx = (uxip1jk - uxijk) / dx;
				variable duydy = (uyijp1k - uyijk) / dy;
				variable duzdz = (uzijkp1 - uzijk) / dz;
				variable tr_eps = duxdx + duydy + duzdz;
				R[Isxx(i, j, k)] = sxx(i, j, k) - ( 2 * mu * duxdx + lam * tr_eps );
				R[Isyy(i, j, k)] = syy(i, j, k) - ( 2 * mu * duydy + lam * tr_eps );
				R[Iszz(i, j, k)] = szz(i, j, k) - ( 2 * mu * duzdz + lam * tr_eps );

				//R[Isxx(i,j,k)] /= E;
				//R[Isyy(i,j,k)] /= E;
				//R[Iszz(i,j,k)] /= E;
			}
		}
	}

	// --------------------------------- Node loop
	// Equations for sxy, sxz, syz:
	//
	// sigma = 2*mu*eps + lam*tr(eps)*I (diagonal parts)
	//
	for (int i = 0; i < N+1; i++) {
		for (int j = 0; j < N+1; j++) {
			for (int k = 0; k < N+1; k++){
				// For dux/dy, duy/dx
				// we need: ux(i,j,k), ux(i,j-1,k), uy(i,j,k), uy(i-1,j,k)
				variable uxijk = 0.0, uxijm1k = 0.0, uyijk = 0.0, uyim1jk = 0.0;
				if(k < N){				
					if (j < N)
						uxijk = ux(i,j,k);
					if (j > 0)
						uxijm1k = ux(i,j-1,k);
					if (i < N)
						uyijk = uy(i,j,k);
					if (i > 0)
						uyim1jk = uy(i-1,j,k);
				}

				// For dux/dz, duz/dx
				// we need: ux(i,j,k), ux(i,j,k-1), uz(i,j,k), uz(i-1,j,k)
				variable uxijkm1 = 0.0, uzijk = 0.0, uzim1jk = 0.0;
				if(j < N){
					if (k > 0)
						uxijkm1 = ux(i,j,k-1);
					if (k < N && i < N) // why??
						uzijk = uz(i,j,k);
					if (i > 0)
						uzim1jk = uz(i-1,j,k);
				}

				// For duy/dz, duz/dy
				// we need: uy(i,j,k), uy(i,j,k-1), uz(i,j,k), uz(i,j-1,k)
				variable uyijkm1 = 0.0, uzijm1k = 0.0;
				if(i < N){
					if (k > 0)
						uyijkm1 = uy(i,j,k-1);
					if (j > 0)
						uzijm1k = uz(i,j-1,k);
				}

				variable duxdy = (uxijk - uxijm1k) / dy;
				variable duydx = (uyijk - uyim1jk) / dx;
				variable duxdz = (uxijk - uxijkm1) / dz;
				variable duzdx = (uzijk - uzim1jk) / dx;
				variable duydz = (uyijk - uyijkm1) / dz;
				variable duzdy = (uzijk - uzijm1k) / dy;
				R[Isxy(i, j, k)] = sxy(i, j, k) - mu * (duxdy + duydx); // sxy = 2*mu*eps_xy == mu * (dux/dy + duy/dx)
				R[Isxz(i, j, k)] = sxz(i, j, k) - mu * (duxdz + duzdx); // sxz = 2*mu*eps_xz == mu * (dux/dz + duz/dx)
				R[Isyz(i, j, k)] = syz(i, j, k) - mu * (duydz + duzdy); // syz = 2*mu*eps_yz == mu * (duy/dz + duz/dy)
				R[Isyx(i, j, k)] = syx(i, j, k) - mu * (duxdy + duydx); // syx = 2*mu*eps_yx == mu * (dux/dy + duy/dx)
				R[Iszx(i, j, k)] = szx(i, j, k) - mu * (duxdz + duzdx); // szx = 2*mu*eps_zx == mu * (dux/dz + duz/dx)
				R[Iszy(i, j, k)] = szy(i, j, k) - mu * (duydz + duzdy); // szy = 2*mu*eps_zy == mu * (duy/dz + duz/dy)
				
				//R[Isxy(i,j,k)] /= E;
				//R[Isxz(i,j,k)] /= E;
				//R[Isyz(i,j,k)] /= E;
			}
		}
	}

	// --------------------------------- Face-x loop
	// Equations for u_x
	//
	// div sigma = -F     (x-coordinate, F_x = 0)
	for (int i = 0; i < N+1; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++){
				if (i > 0 && i < N && j > 0 && j < N-1 && k > 0 && k < N-1) {
					variable dsxxdx = (sxx(i, j,   k  ) - sxx(i-1, j, k)) / dx;
					variable dsyxdy = (syx(i, j+1, k  ) - syx(i,   j, k)) / dy;
					variable dszxdz = (szx(i, j,   k+1) - szx(i,   j, k)) / dz;
					R[Iux(i, j, k)] = dsxxdx + dsyxdy + dszxdz - 0.0;
				}
				else
					R[Iux(i, j, k)] = ux(i,j,k);
				R[Iux(i, j, k)] *= -1;
			}
				
		}
	}

	// --------------------------------- Face-y loop
	// Equations for u_y
	//
	// div sigma = -F     (y-coordinate, F_y = 0)
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N+1; j++) {
			for (int k = 0; k < N; k++){
				if (i > 0 && i < N-1 && j > 0 && j < N && k > 0 && k < N-1) {
					variable dsxydx = (sxy(i+1, j, k  ) - sxy(i, j,   k)) / dx;
					variable dsyydy = (syy(i,   j, k  ) - syy(i, j-1, k)) / dy;
					variable dszydz = (szy(i,   j, k+1) - szy(i, j,   k)) / dz;
					R[Iuy(i, j, k)] = dsxydx + dsyydy + dszydz - 0.0;
				}
				else
					R[Iuy(i, j, k)] = uy(i,j,k);
				R[Iuy(i, j, k)] *= -1;
				}
		}
	}

	// --------------------------------- Face-z loop
	// Equations for u_z
	//
	// div sigma = -F     (z-coordinate, F_z = -1)
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N+1; k++){
				if (i > 0 && i < N-1 && j > 0 && j < N-1 && k > 0 && k < N) {
					variable dsxzdx = (sxz(i+1, j,   k) - sxz(i, j, k  )) / dx;
					variable dsyzdy = (syz(i,   j+1, k) - syz(i, j, k  )) / dy;
					variable dszzdz = (szz(i,   j,   k) - szz(i, j, k-1)) / dz;
					R[Iuz(i, j, k)] = dsxzdx + dsyzdy + dszzdz - 1.0;
				}
				else
					R[Iuz(i, j, k)] = uz(i,j,k);
				R[Iuz(i, j, k)] *= -1;
				}
		}
	}
	//std::cout << "System is assembled" << std::endl;
}

void Problem::run()
{
	// Total number of unknowns:
	// Cell   ( N   * N   ):     2
	// Node   ((N+1)*(N+1)):     1
	// Face_x ( N   *(N+1)):     1
	// Face_y ( N   *(N+1)):     1
	int tot_size = N*N*N*3 + N*N*(N+1)*3 + (N+1)*(N+1)*(N+1)*6;
	std::cout << "Total number of unknowns: " << tot_size << std::endl;
	Residual R("residual", 0, tot_size);
	sol = Sparse::Vector("solution", 0, tot_size);
	Sparse::Vector update = Sparse::Vector("newton_update", 0, tot_size);

	//Solver S("trilinos_aztec");
	//S.SetParameter("drop_tolerance", "2e0");
	Solver S("inner_mptiluc");
	S.SetParameter("drop_tolerance", "1e-1");//std::to_string(N/80));
	S.SetParameter("reuse_tolerance", "1e3");
	S.SetParameter("absolute_tolerance", "1e-15");
	S.SetParameter("relative_tolerance", "1e-10");
	S.SetParameter("maximum_iterations", "10000");
	S.SetParameter("verbosity", "3");


	// ============================== Newton loop
	std::cout << std::endl << "Starting Newton loop" << std::endl;
	double r = 1.0, r0 = 1.0;
	bool converged = false;
	int maxit = 20;
	double rtol = 1e-6, atol = 1e-8, divtol = 1e10;
	for(int nit = 0; nit < maxit; nit++){
		fillResidual(R);
		//R.GetJacobian().Save("J.mtx");

		// Convergence check
		r = R.Norm();
		if(nit == 0)
			r0 = r;
		std:: cout << "  iter " << nit << ", |r|_2 = " << r << std::endl;
		if(r < atol || r < rtol * r0){
			std::cout << "Newton converged!" << std:: endl;
			converged = true;
			break;
		}
		if(r > divtol){
			std::cout << "Newton diverged!" << std:: endl;
			break;
		}

		S.SetMatrix(R.GetJacobian());
		bool solved = S.Solve(R.GetResidual(), update);
		if (!solved){
			std::cout << "Linear solver failed: " << S.GetReason() << std::endl;
			std::cout << "Residual: " << S.Residual() << std::endl;
			exit(-1);
		}
		//std::cout << "Lin.it:   " << S.Iterations() << std::endl;
		//std::cout << "Residual: " << S.Residual() << std::endl;

		//double solmax = 0.0;
		for (unsigned i = 0; i < sol.Size(); i++) {
			sol[i] -= update[i];
			//solmax = std::max(solmax, abs(sol[i]));
		}
		//std::cout << "Max. abs. val. in sol = " << solmax << std::endl;
	}
	if(!converged){
		std::cout << "Newton failed to converge" << std::endl;
		exit(-1);
	}
}

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	if (argc < 2) {
		std::cout << "Usage: mech2d <N>" << std::endl;
		exit(-1);
	}

	Problem P(atoi(argv[1]));
	P.run();
	//P.saveVTK2D();
	P.saveVTK3D();

	std::cout << "Success!" << std::endl;
	return 0;
}
