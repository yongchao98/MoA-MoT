import numpy as np
from scipy.optimize import minimize

def solve_sphere_packing(N, R_cyl, H_cyl, r_sph):
    """
    Attempts to find a feasible packing for N spheres in a cylinder.

    This function formulates the sphere packing problem as a numerical
    optimization problem and uses SciPy's SLSQP solver to find a
    set of sphere coordinates that satisfy the geometric constraints.

    Args:
        N (int): The number of spheres to pack.
        R_cyl (float): The radius of the cylinder.
        H_cyl (float): The height of the cylinder.
        r_sph (float): The radius of the spheres.
    """
    print(f"--- Attempting to pack {N} spheres ---")
    print("Problem Parameters:")
    print(f"  Cylinder: Radius={R_cyl}, Height={H_cyl}")
    print(f"  Sphere: Radius={r_sph}")
    print("----------------------------------------")

    # Effective radius and height for sphere centers
    R_eff = R_cyl - r_sph
    Z_min = r_sph
    Z_max = H_cyl - r_sph
    min_dist_sq = (2 * r_sph)**2

    # The solver operates on a single flat array of variables p = [x1,y1,z1,x2,y2,z2,...]
    # Objective function: We don't need to optimize anything, just find a
    # feasible point, so the objective is constant zero.
    def objective(p):
        return 0

    # --- Define Constraints ---
    constraints = []

    # 1. Non-overlap constraints for every pair of spheres
    # (xi-xj)^2 + (yi-yj)^2 + (zi-zj)^2 - (2r)^2 >= 0
    for i in range(N):
        for j in range(i + 1, N):
            # The lambda function captures the current i and j
            con_fun = lambda p, i=i, j=j: \
                (p[i*3] - p[j*3])**2 + \
                (p[i*3+1] - p[j*3+1])**2 + \
                (p[i*3+2] - p[j*3+2])**2 - min_dist_sq
            constraints.append({'type': 'ineq', 'fun': con_fun})

    # 2. Cylinder containment constraints for each sphere
    # (R_eff)^2 - (xi^2 + yi^2) >= 0
    for i in range(N):
        con_fun = lambda p, i=i: R_eff**2 - (p[i*3]**2 + p[i*3+1]**2)
        constraints.append({'type': 'ineq', 'fun': con_fun})
    
    # --- Define Bounds for each variable ---
    # xi, yi in [-R_eff, R_eff]
    # zi in [Z_min, Z_max]
    bounds = []
    for i in range(N):
        bounds.extend([(-R_eff, R_eff), (-R_eff, R_eff), (Z_min, Z_max)])

    # --- Create a good initial guess (p0) ---
    # Randomly distribute points within the allowed volume for sphere centers.
    # This is better than a simple grid, which can start in an invalid state.
    p0 = np.random.rand(N * 3)
    p0_reshaped = p0.reshape((N, 3))
    
    # z-coordinates
    p0_reshaped[:, 2] = p0_reshaped[:, 2] * (Z_max - Z_min) + Z_min
    # x,y-coordinates (uniformly in a circle)
    r_rand = np.sqrt(np.random.rand(N)) * R_eff
    theta_rand = np.random.rand(N) * 2 * np.pi
    p0_reshaped[:, 0] = r_rand * np.cos(theta_rand)
    p0_reshaped[:, 1] = r_rand * np.sin(theta_rand)
    p0_flat = p0_reshaped.flatten()

    # --- Run the Optimizer ---
    print("Starting optimizer... (This may take a minute)")
    result = minimize(objective, p0_flat, method='SLSQP', bounds=bounds,
                      constraints=constraints, options={'disp': True, 'maxiter': 1000})

    # --- Print Results ---
    print("\n--- Optimization Finished ---")
    if result.success:
        print(f"SUCCESS: A feasible packing for {N} spheres was found.")
        print("Final sphere center coordinates (x, y, z):")
        final_coords = result.x.reshape(N, 3)
        for i in range(N):
            print(f"  Sphere {i+1:2d}: ({final_coords[i, 0]:8.4f}, {final_coords[i, 1]:8.4f}, {final_coords[i, 2]:8.4f})")
    else:
        print(f"FAILURE: Could not find a feasible packing for {N} spheres.")
        print(f"Solver message: {result.message}")

    return result.success

# --- Main Execution ---
if __name__ == "__main__":
    # Problem parameters
    NUM_SPHERES = 32
    CYLINDER_RADIUS = 10.0
    CYLINDER_HEIGHT = 11.0
    SPHERE_RADIUS = 2.0
    
    solve_sphere_packing(NUM_SPHERES, CYLINDER_RADIUS, CYLINDER_HEIGHT, SPHERE_RADIUS)

