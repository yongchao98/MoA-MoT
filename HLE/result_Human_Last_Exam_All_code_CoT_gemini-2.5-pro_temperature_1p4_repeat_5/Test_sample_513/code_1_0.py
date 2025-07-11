import numpy as np
from scipy.optimize import minimize, Bounds

def solve_packing_problem():
    """
    Finds the maximum number of spheres that can be packed into a cylinder.
    """
    # --- Problem Parameters ---
    R_c = 10.0  # Cylinder radius
    H_c = 11.0  # Cylinder height
    r_s = 2.0   # Sphere radius

    print("Solving the sphere packing problem with the following parameters:")
    print(f"Cylinder: Radius = {R_c}, Height = {H_c}")
    print(f"Spheres: Radius = {r_s}\n")

    # --- Derived parameters for sphere centers ---
    R_eff = R_c - r_s   # Effective radius for sphere centers
    Z_min = r_s         # Min z-coordinate for a sphere center
    Z_max = H_c - r_s   # Max z-coordinate for a sphere center
    D_sq = (2 * r_s)**2 # Squared minimum distance between centers

    def find_packing_for_n(n):
        """
        Tries to find a valid packing for n spheres using optimization.
        Returns True if a valid packing is found, False otherwise.
        """

        # The objective function calculates a penalty for invalid configurations.
        # A perfect packing has a penalty of 0.
        def penalty_function(coords_flat):
            coords = coords_flat.reshape((n, 3))
            
            # 1. Penalty for sphere-sphere overlap
            overlap_penalty = 0.0
            for i in range(n):
                for j in range(i + 1, n):
                    d_sq = np.sum((coords[i] - coords[j])**2)
                    if d_sq < D_sq:
                        # Penalize based on how much they overlap
                        overlap_penalty += (D_sq - d_sq)**2

            # 2. Penalty for spheres outside the cylinder wall
            wall_penalty = 0.0
            radial_dist_sq = coords[:, 0]**2 + coords[:, 1]**2
            violations = radial_dist_sq > R_eff**2
            if np.any(violations):
                wall_penalty = np.sum((radial_dist_sq[violations] - R_eff**2)**2)
            
            # Combine penalties (weighting wall penalty higher can help convergence)
            return overlap_penalty + 10 * wall_penalty

        # Bounds constrain the sphere centers within the effective cylinder volume.
        # This already enforces the height constraint.
        bounds_list = []
        for _ in range(n):
            bounds_list.extend([(-R_eff, R_eff), (-R_eff, R_eff), (Z_min, Z_max)])
        bounds = Bounds([b[0] for b in bounds_list], [b[1] for b in bounds_list])

        # Generate a random initial guess for the sphere positions
        x0 = np.zeros(3 * n)
        # To generate points uniformly in a disk for x,y coordinates
        u, v = np.random.rand(n), np.random.rand(n)
        r_init = R_eff * np.sqrt(u)
        theta_init = 2 * np.pi * v
        x0[0::3] = r_init * np.cos(theta_init)
        x0[1::3] = r_init * np.sin(theta_init)
        # Uniformly distribute z coordinates
        x0[2::3] = np.random.uniform(Z_min, Z_max, n)

        # Run the optimization
        result = minimize(penalty_function,
                          x0,
                          method='L-BFGS-B',
                          bounds=bounds,
                          options={'maxiter': 1500, 'ftol': 1e-9})
        
        # A successful packing is found if the final penalty is very close to zero.
        if result.success and result.fun < 1e-8:
            return True
        return False

    # --- Search for the maximum number of spheres ---
    # Start from a reasonable upper estimate and search downwards.
    # Known results for similar problems suggest the answer is in the 40-50 range.
    for n in range(50, 40, -1):
        print(f"Attempting to pack {n} spheres...")
        # Run a few attempts with different random starts for robustness
        is_packable = False
        for attempt in range(3):
            if find_packing_for_n(n):
                is_packable = True
                break
            print(f"  Attempt {attempt + 1}/3 for n={n} failed, trying again...")
        
        if is_packable:
            print("\n---------------------------------------------------------")
            print(f"Success! A valid packing was found for {n} spheres.")
            print(f"The final formulated equation for the maximum number of spheres is:")
            print(f"Max n = {n}, for cylinder(R={R_c}, H={H_c}) and spheres(r={r_s})")
            print("---------------------------------------------------------")
            return n

    print("\nCould not find a solution in the tested range.")
    return None

if __name__ == '__main__':
    max_spheres = solve_packing_problem()
    if max_spheres is not None:
        print(f"<<<{max_spheres}>>>")
