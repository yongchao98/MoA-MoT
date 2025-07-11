import numpy as np
from scipy.optimize import minimize
from scipy.spatial.distance import pdist

def solve_packing_problem():
    """
    Solves the sphere packing problem in a cylinder using optimization.
    """
    # Problem parameters
    R_cyl = 10.0
    H_cyl = 11.0
    r_sph = 2.0

    # Derived parameters for the centers' allowed positions
    R_eff = R_cyl - r_sph
    z_min = r_sph
    z_max = H_cyl - r_sph
    d_sq = (2 * r_sph)**2  # Squared diameter for distance checks

    def check_feasibility(n_spheres, n_trials=10, max_iter=2000):
        """
        Tries to find a feasible packing for n_spheres by minimizing a penalty function.
        Returns True if a valid packing is found, False otherwise.
        """
        print(f"Checking feasibility for {n_spheres} spheres...")

        def objective_function(p):
            """
            Penalty function. Its value is 0 for a valid packing, >0 otherwise.
            p is a flat array of coordinates (x1, y1, z1, x2, y2, z2, ...).
            """
            centers = p.reshape((n_spheres, 3))

            # Penalty 1: Overlap between spheres
            # pdist calculates the pairwise distances between all points
            sq_dists = pdist(centers, 'sqeuclidean')
            overlap_violations = np.maximum(0, d_sq - sq_dists)
            penalty_overlap = np.sum(overlap_violations)

            # Penalty 2: Spheres outside the cylinder's radial boundary
            radial_sq_pos = centers[:, 0]**2 + centers[:, 1]**2
            radial_violations = np.maximum(0, radial_sq_pos - R_eff**2)
            penalty_radial = np.sum(radial_violations)
            
            return penalty_overlap + penalty_radial

        # Bounds for each variable (x_i, y_i, z_i)
        bounds = []
        for _ in range(n_spheres):
            bounds.extend([(-R_eff, R_eff), (-R_eff, R_eff), (z_min, z_max)])

        # Run optimization multiple times from random starting points
        for i in range(n_trials):
            # Generate a random initial guess for the sphere centers
            p0 = np.random.rand(n_spheres * 3)
            p0_reshaped = p0.reshape((n_spheres, 3))
            
            # Scale initial guess to fit within the bounds
            p0_reshaped[:, 0] = p0_reshaped[:, 0] * 2 * R_eff - R_eff
            p0_reshaped[:, 1] = p0_reshaped[:, 1] * 2 * R_eff - R_eff
            p0_reshaped[:, 2] = p0_reshaped[:, 2] * (z_max - z_min) + z_min
            
            # Run the L-BFGS-B optimizer
            res = minimize(
                objective_function,
                p0,
                method='L-BFGS-B',
                bounds=bounds,
                options={'ftol': 1e-9, 'maxiter': max_iter}
            )

            # If penalty is virtually zero, a solution is found
            if res.success and res.fun < 1e-6:
                print(f"  SUCCESS: Feasible packing found for {n_spheres} spheres (penalty={res.fun:.2e}).")
                return True
        
        print(f"  FAILURE: No feasible packing found for {n_spheres} spheres after {n_trials} trials.")
        return False

    # Binary search for the maximum number of spheres N
    # Lower bound from simple packing estimates, upper bound from volume considerations.
    low = 24  # A simple 2-layer packing
    high = 50 # A safe upper bound
    best_n = 0

    while low <= high:
        mid = (low + high) // 2
        if mid == 0:
            break
        if check_feasibility(mid):
            best_n = mid  # This is a possible solution
            low = mid + 1 # Try for more
        else:
            high = mid - 1 # Too many, try for less
    
    print("\n" + "="*40)
    print("Final Problem Specification and Result")
    print("="*40)
    print(f"Sphere Radius (r): {r_sph}")
    print(f"Cylinder Radius (R): {R_cyl}")
    print(f"Cylinder Height (H): {H_cyl}")
    print("-" * 40)
    print(f"The maximum number of spheres that can be packed is: {best_n}")
    print("="*40)

if __name__ == '__main__':
    solve_packing_problem()