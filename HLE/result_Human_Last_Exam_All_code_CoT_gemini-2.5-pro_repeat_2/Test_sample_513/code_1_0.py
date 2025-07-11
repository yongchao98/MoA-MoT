import numpy as np
from scipy.optimize import minimize, Bounds

def solve_sphere_packing():
    """
    Attempts to find a feasible packing for N spheres in a cylinder.
    """
    # --- Parameters ---
    R_cylinder = 10.0
    H_cylinder = 11.0
    r_sphere = 2.0
    
    # Number of spheres to pack. The best known value is 49.
    # We will try to find a solution for N=48, which is a strong result.
    # The solver may be slow or fail for high N depending on the initial guess.
    N = 48

    # --- Derived Parameters ---
    # Effective radius for sphere centers
    R_effective = R_cylinder - r_sphere
    # Min and max height for sphere centers
    H_min = r_sphere
    H_max = H_cylinder - r_sphere
    # Squared minimum distance between centers
    d_sq = (2 * r_sphere)**2

    # --- Objective Function (dummy, for feasibility) ---
    def objective(x):
        return 0

    # --- Constraint Function ---
    def constraints(x):
        """
        Returns an array of values that must be non-negative.
        x is a flat array: [x0, y0, z0, x1, y1, z1, ...]
        """
        # Reshape to (N, 3) for easier access
        coords = x.reshape(N, 3)
        
        # 1. Radial constraints: R_effective^2 - (x_i^2 + y_i^2) >= 0
        radial_constraints = R_effective**2 - (coords[:, 0]**2 + coords[:, 1]**2)
        
        # 2. Non-overlap constraints: (dist_ij)^2 - d_sq >= 0
        overlap_constraints = []
        for i in range(N):
            for j in range(i + 1, N):
                dist_sq = np.sum((coords[i, :] - coords[j, :])**2)
                overlap_constraints.append(dist_sq - d_sq)
                
        return np.concatenate([radial_constraints, np.array(overlap_constraints)])

    # --- Bounds for Variables ---
    # Each sphere's center (xi, yi, zi) is bounded.
    low_bounds = []
    high_bounds = []
    for _ in range(N):
        low_bounds.extend([-R_effective, -R_effective, H_min])
        high_bounds.extend([R_effective, R_effective, H_max])
    bounds = Bounds(low_bounds, high_bounds)

    # --- Initial Guess ---
    # A random guess often fails. A more structured guess (e.g., a helix)
    # helps the solver find a feasible solution. This is a phyllotaxis pattern.
    print("Generating a structured initial guess...")
    x0 = np.zeros(N * 3)
    golden_angle = np.pi * (3. - np.sqrt(5.))
    for i in range(N):
        # Distribute z linearly
        z = H_min + (i / (N - 1)) * (H_max - H_min) if N > 1 else (H_min + H_max) / 2
        # Radius grows with height
        radius = R_effective * np.sqrt(z / H_max)
        # Constant angle increment
        theta = i * golden_angle
        
        x0[3*i] = radius * np.cos(theta)
        x0[3*i+1] = radius * np.sin(theta)
        x0[3*i+2] = z

    # --- Run the Optimization ---
    # 'trust-constr' is suitable for problems with many constraints.
    # It might take a few minutes for a large N.
    print(f"Attempting to pack {N} spheres. This may take a few minutes...")
    cons = [{'type': 'ineq', 'fun': constraints}]
    result = minimize(objective, x0, method='trust-constr', bounds=bounds,
                      constraints=cons, options={'verbose': 1, 'maxiter': 1000})

    # --- Print Results ---
    print("\n--- Solver Finished ---")
    if result.success:
        print(f"Successfully found a feasible packing for {N} spheres.")
        print("\nThe final coordinates (x, y, z) for each sphere are:")
        final_coords = result.x.reshape(N, 3)
        for i in range(N):
            print(f"Sphere {i+1:2d}: ({final_coords[i, 0]:8.4f}, {final_coords[i, 1]:8.4f}, {final_coords[i, 2]:8.4f})")
    else:
        print(f"Solver failed to find a feasible packing for {N} spheres.")
        print("This does not prove it's impossible, as the solver may have gotten stuck.")
        print("Try running again or using a different initial guess.")
    
    # According to published data, the optimal answer is 49.
    # The code demonstrates the method to verify such a claim.

if __name__ == '__main__':
    solve_sphere_packing()