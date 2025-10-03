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
//   Cell   index range: 0...N-1, 0...N-1
//   Node   index range: 0...N,   0...N
//   Face-x index range: 0...N,   0...N-1
//   Face-y index range: 0...N-1, 0...N
//
//   Cell   vars: sigma_xx, sigma_yy
//   Node   vars: sigma_xy
//   Face-x vars: u_x
//   Face-y vars: u_y

class Problem
{
private:
	// =====================================================================================================
	// Problem parameters: discretization and physical parameters
	int N;
	double dx;
	double dy;
	const double Lx = 1.0;
	const double Ly = Lx;
	const double E = 3.5e6;                               // Young's modulus, MPa
	const double nu = 0.3;                                // Poisson ratio
	const double lam = E * nu / (1 + nu) / (1 - 2 * nu);  // Lame parameter lambda
	const double mu = E / 2 / (1 + nu);                   // Lame parameter mu


	// =====================================================================================================
	// INMOST functionality
	/// Solution vector that contains all the unknowns: sigma_xx, sigma_yy, sigma_xy, u_x, u_y
	INMOST::Sparse::Vector sol;
	/// Linear solver for systems with Jacobian matrix
	//INMOST::Solver S;
	/// Nonlinear residual
	//INMOST::Residual R;


	// =====================================================================================================
	// Index functions that map (i,j) for a specific unknown to the position in the global vector 'sol'
	int Isxx(int i, int j); // sigma_xx
	int Isyy(int i, int j); // sigma_yy
	int Isxy(int i, int j); // sigma_xy
	int Iux(int i, int j);  // u_x
	int Iuy(int i, int j);  // u_y


	// =====================================================================================================
	// Functions that construct unknonwns from given locations in 'sol'
	unknown sxx(int i, int j) { return unknown(sol[Isxx(i,j)], Isxx(i,j)); }
	unknown syy(int i, int j) { return unknown(sol[Isyy(i,j)], Isyy(i,j)); }
	unknown sxy(int i, int j) { return unknown(sol[Isxy(i,j)], Isxy(i,j)); }
	unknown ux(int i, int j)  { return unknown(sol[Iux(i,j)],  Iux(i,j));  }
	unknown uy(int i, int j)  { return unknown(sol[Iuy(i,j)],  Iuy(i,j));  }



public:
	Problem(int N_)
	{
		N = N_;
		dx = Lx / N;
		dy = Ly / N;
	}
	~Problem() {}
	void run();
	void fillResidual(Residual &R);
	void saveVTK();
};


int Problem::Isxx(int i, int j)
{
	if (i < 0 || i > N-1) 
		std::cout << "Isxx: wrong i = " << i << " for N = " << N << std::endl;
	if (j < 0 || j > N-1)
		std::cout << "Isxx: wrong j = " << j << " for N = " << N << std::endl;
	return i*N + j;
}

int Problem::Isyy(int i, int j)
{
	if (i < 0 || i > N-1)
		std::cout << "Isyy: wrong i = " << i << " for N = " << N << std::endl;
	if (j < 0 || j > N-1)
		std::cout << "Isyy: wrong j = " << j << " for N = " << N << std::endl;
	return i*N + j + N*N;
}

int Problem::Isxy(int i, int j)
{
	if (i < 0 || i > N)
		std::cout << "Isxy: wrong i = " << i << " for N = " << N << std::endl;
	if (j < 0 || j > N)
		std::cout << "Isxy: wrong j = " << j << " for N = " << N << std::endl;
	return i*(N+1) + j + N*N*2;
}

int Problem::Iux(int i, int j)
{
	if (i < 0 || i > N)
		std::cout << "Iux: wrong i = " << i << " for N = " << N << std::endl;
	if (j < 0 || j > N-1)
		std::cout << "Iux: wrong j = " << j << " for N = " << N << std::endl;
	return (i)*N + j + N * N * 2 + (N + 1) * (N + 1);
}

int Problem::Iuy(int i, int j)
{
	if (i < 0 || i > N-1)
		std::cout << "Iuy: wrong i = " << i << " for N = " << N << std::endl;
	if (j < 0 || j > N)
		std::cout << "Iuy: wrong j = " << j << " for N = " << N << std::endl;
	return ((i) * (N + 1) + j + N * N * 2 + (N + 1) * (N + 1) + N * (N + 1));
}

void Problem::saveVTK()
{
	
	// ===================================================================================
	// Save VTK with results (legacy code)
	std::ofstream out;
	out.open("res.vtk");
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
				uxij = sol[Iux(i, j)];
			if (i < nx && j > 0)
				uxijm1 = sol[Iux(i, j - 1)];
			if (i < nx && j < ny)
				uyij = sol[Iuy(i, j)];
			if (i > 0 && j < ny)
				uyim1j = sol[Iuy(i - 1, j)];

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
			out << sol[Isyy(i,j)] << std::endl;
		}
	}

	out << "SCALARS Stress_xx double" << std::endl;
	out << "LOOKUP_TABLE default" << std::endl;
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			out << sol[Isxx(i,j)] << std::endl;
		}
	}

	// Values in cell are averaged from cell nodes
	out << "SCALARS Txy double" << std::endl;
	out << "LOOKUP_TABLE default" << std::endl;
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			//out << 0.25 * (Txy[i + j * (nx + 1)] + Txy[i + j * (nx + 1)]
			//	+ Txy[i + 1 + j * (nx + 1)] + Txy[i + (j + 1) * (nx + 1)]) << std::endl;
			out << 0.25 * (sol[Isxy(i,j)] + sol[Isxy(i+1, j)]
				+ sol[Isxy(i, j+1)] + sol[Isxy(i+1, j+1)]) << std::endl;
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
				uxij = sol[Iux(i,j)]; 
			if (i < nx && j > 0)
				uxijm1 = sol[Iux(i,j-1)];
			if (i < nx && j < ny)
				uyij = sol[Iuy(i,j)];
			if (i > 0 && j < ny)
				uyim1j = sol[Iuy(i-1,j)];

			double ux = 0.5 * (uxij + uxijm1);
			double uy = 0.5 * (uyij + uyim1j);

			out <<  ux << " " << uy << " 0.0" << std::endl;
		}
	}

	out.close();
}

void Problem::fillResidual(Residual &R)
{
	R.Clear();
	// --------------------------------- Cell loop
	// Equations for sxx, syy:
	// 
	// sigma = C u (diagonal parts)
	//
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			variable uxij = 0.0, uxip1j = 0.0, uyij = 0.0, uyijp1 = 0.0;
			if (i > 0)
				uxij = ux(i, j);
			if (i < N+1)
				uxip1j = ux(i+1, j);
			if (j > 0)
				uyij = uy(i, j);
			if (j < N+1)
				uyijp1 = uy(i, j+1);
			variable duxdx = (uxip1j - uxij) / dx;
			variable duydy = (uyijp1 - uyij) / dy;
			R[Isxx(i, j)] = sxx(i, j) - (2 * mu + lam) * duxdx + lam * duydy;
			R[Isyy(i, j)] = syy(i, j) - (2 * mu + lam) * duydy + lam * duxdx;
		}
	}

	// --------------------------------- Node loop
	// Equations for sxy:
	//
	// sigma = C u (diagonal parts)
	//
	for (int i = 0; i < N + 1; i++) {
		for (int j = 0; j < N+1; j++) {
			// For dux/dy, duy/dx
			// we need: ux(i,j), ux(i,j-1), uy(i,j), uy(i-1,j)
			variable uxij = 0.0, uxijm1 = 0.0, uyij = 0.0, uyim1j = 0.0;
			if (j < N)
				uxij = ux(i,j);
			if (j > 0)
				uxijm1 = ux(i,j-1);
			if (i < N)
				uyij = uy(i,j);
			if (i > 0)
				uyim1j = uy(i-1,j);

			variable duxdy = (uxij - uxijm1) / dy;
			variable duydx = (uyij - uyim1j) / dx;
			R[Isxy(i, j)] = sxy(i, j) - mu * (duxdy + duydx);
			//R[Isxy(i, j)] = sxy(i, j) - 3;
		}
	}

	// --------------------------------- Face-x loop
	// Equations for u_x
	//
	// div sigma = -F     (x-coordinate, F_x = 0)
	for (int i = 0; i < N+1; i++) {
		for (int j = 0; j < N; j++) {
			if (i > 0 && i < N && j > 0 && j < N-1) {
				variable dsxxdx = (sxx(i, j  ) - sxx(i-1, j)) / dx;
				variable dsxydy = (sxy(i, j+1) - sxy(i,   j)) / dy;
				R[Iux(i, j)] = dsxxdx + dsxydy - 0.0;
			}
			else
				R[Iux(i, j)] = ux(i,j);
		}
	}

	// --------------------------------- Face-y loop
	// Equations for u_y
	//
	// div sigma = -F     (y-coordinate, F_y = -1)
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N+1; j++) {
			if (j > 0 && j < N && i > 0 && i < N-1) {
				variable dsyydy = (syy(i,   j) - syy(i, j-1)) / dy;
				variable dsxydx = (sxy(i+1, j) - sxy(i, j  )) / dx;
				R[Iuy(i, j)] = dsyydy + dsxydx - 1.0;
				double x = (i + 0.5) * dx;
				double y = j * dy;
				//R[Iuy(i, j)] = uy(i, j);// -x * (1 - x) * y * (1 - y);
			}
			else
				R[Iuy(i, j)] = uy(i,j);
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
	int tot_size = N*N*2 + N*(N+1)*2 + (N+1)*(N+1);
	std::cout << "Total number of unknowns: " << tot_size << std::endl;
	Residual R("residual", 0, tot_size);
	sol = Sparse::Vector("solution", 0, tot_size);
	Sparse::Vector update = Sparse::Vector("newton_update", 0, tot_size);

	Solver S("inner_mptiluc");
	S.SetParameter("drop_tolerance", "1e-6");
	S.SetParameter("absolute_tolerance", "1e-15");
	S.SetParameter("relative_tolerance", "1e-10");
	S.SetParameter("maximum_iterations", "10000");


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
	P.saveVTK();

	std::cout << "Success!" << std::endl;
	return 0;
}