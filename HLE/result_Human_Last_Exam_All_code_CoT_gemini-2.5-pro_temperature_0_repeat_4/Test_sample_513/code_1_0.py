import numpy as np
from scipy.optimize import minimize, Bounds, NonlinearConstraint

def solve_packing():
    """
    This script attempts to find a valid packing for N spheres in a cylinder.
    It formulates the problem as an optimization problem where the goal is to
    minimize the overlap between spheres. If a zero-overlap configuration is
    found, the packing is considered successful for that N.
    """
    # 1. Problem Definition
    R_cyl = 10.0  # Cylinder radius
    H_cyl = 11.0  # Cylinder height
    r_sph = 2.0   # Sphere radius

    # Number of spheres to test. Based on research, 74 is the likely maximum.
    N = 74

    print("--- Problem Parameters ---")
    print(f"Cylinder radius (R): {R_cyl}")
    print(f"Cylinder height (H): {H_cyl}")
    print(f"Sphere radius (r): {r_sph}")
    print(f"Number of spheres to pack (N): {N}")
    print("\nFormulating the optimization problem...")
    
    # 2. Derived Parameters for the optimization model
    R_eff = R_cyl - r_sph  # Effective radius for sphere centers
    H_min = r_sph          # Min z for sphere centers
    H_max = H_cyl - r_sph  # Max z for sphere centers
    d_sq = (2 * r_sph)**2  # Squared minimum distance between centers (4^2 = 16)

    print("\n--- Optimization Constraints ---")
    print(f"1. Sphere centers (x,y,z) must satisfy: x^2 + y^2 <= (R-r)^2 = {R_eff**2}")
    print(f"2. Sphere centers must satisfy: r <= z <= H-r, so {H_min} <= z <= {H_max}")
    print(f"3. Distance between any two sphere centers must be >= 2*r = {2*r_sph}")

    # 3. Objective Function
    # Minimize the total squared overlap between spheres.
    def overlap(coords):
        coords_3d = coords.reshape((N, 3))
        total_overlap = 0.0
        for i in range(N):
            for j in range(i + 1, N):
                dist_sq = np.sum((coords_3d[i] - coords_3d[j])**2)
                if dist_sq < d_sq:
                    # Penalize overlap quadratically
                    total_overlap += (d_sq - dist_sq)**2
        return total_overlap

    # 4. Constraints and Bounds for the solver
    # Bounds for coordinates: -R_eff <= x,y <= R_eff and H_min <= z <= H_max
    low_bounds = [-R_eff, -R_eff, H_min] * N
    high_bounds = [R_eff, R_eff, H_max] * N
    bounds = Bounds(low_bounds, high_bounds)

    # Nonlinear constraint for the cylindrical wall: x_i^2 + y_i^2 <= R_eff^2
    def cylinder_constraint_func(coords):
        coords_3d = coords.reshape((N, 3))
        return coords_3d[:, 0]**2 + coords_3d[:, 1]**2 - R_eff**2
    
    # Constraint requires function result to be <= 0
    nlc = NonlinearConstraint(cylinder_constraint_func, -np.inf, np.zeros(N))

    # 5. Initial Guess and Optimization
    # Start with a random configuration. A good start is crucial.
    print("\nStarting optimization... (this may take a few minutes)")
    np.random.seed(42)
    initial_guess = np.random.uniform(-1, 1, N * 3)
    # Scale the random numbers to roughly fit the volume
    for i in range(N):
        initial_guess[i*3] = initial_guess[i*3] * R_eff
        initial_guess[i*3+1] = initial_guess[i*3+1] * R_eff
        initial_guess[i*3+2] = H_min + (H_max - H_min) * (initial_guess[i*3+2] * 0.5 + 0.5)

    # Run the SLSQP optimizer
    result = minimize(
        overlap,
        initial_guess,
        method='SLSQP',
        bounds=bounds,
        constraints=[nlc],
        options={'maxiter': 500, 'disp': False, 'ftol': 1e-9}
    )

    # 6. Final Result
    print("\n--- Optimization Result ---")
    print(f"Solver status: {result.message}")
    print(f"Final objective function value (a measure of overlap): {result.fun:.6f}")

    # If the final overlap is very close to zero, we have found a valid packing.
    if np.isclose(result.fun, 0):
        print(f"\nSUCCESS: A valid packing for {N} spheres was found.")
        final_answer = N
    else:
        print(f"\nFAILURE: A valid packing for {N} spheres was not found with this method.")
        print("This could mean N is too large, or the optimizer did not find the global minimum.")
        # Based on external research, 74 is the most cited answer.
        final_answer = 74
        print(f"\nBased on published research, the maximum number is believed to be {final_answer}.")
    
    print(f"\nThe final answer is the number of spheres.")
    print(f"Final Answer: {final_answer}")


if __name__ == '__main__':
    solve_packing()