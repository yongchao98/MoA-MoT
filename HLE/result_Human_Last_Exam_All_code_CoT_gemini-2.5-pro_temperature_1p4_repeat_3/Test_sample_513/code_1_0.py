import numpy as np
from scipy.optimize import minimize

def solve_sphere_packing(num_spheres):
    """
    Formulates and attempts to solve the sphere packing feasibility problem.

    Args:
        num_spheres (int): The number of spheres to try and pack.

    Returns:
        The result from the optimization solver.
    """
    # --- Problem Parameters ---
    R_c = 10.0  # Cylinder radius
    H_c = 11.0  # Cylinder height
    r_s = 2.0   # Sphere radius
    N = num_spheres

    print(f"--- Sphere Packing Problem ---")
    print(f"Objective: Find a feasible packing for N = {N} spheres.")
    print("\nParameters:")
    print(f"  Sphere Radius (r_s): {r_s}")
    print(f"  Cylinder Radius (R_c): {R_c}")
    print(f"  Cylinder Height (H_c): {H_c}")

    # --- Formulation Equations ---
    print("\nFormulation (for any two spheres i and j):")
    # Non-overlap: (xi-xj)^2 + (yi-yj)^2 + (zi-zj)^2 >= (2*rs)^2
    print(f"  1. Non-overlap: (x_i-x_j)^2 + (y_i-y_j)^2 + (z_i-z_j)^2 >= ({2*r_s})^2 = { (2*r_s)**2 }")
    # Radial containment: xi^2 + yi^2 <= (Rc-rs)^2
    print(f"  2. Radial containment: x_i^2 + y_i^2 <= ({R_c}-{r_s})^2 = { (R_c-r_s)**2 }")
    # Height containment: rs <= zi <= Hc-rs
    print(f"  3. Height containment: {r_s} <= z_i <= {H_c-r_s} = {H_c - r_s}")
    print("-" * 30)

    # The problem is to find a point in a high-dimensional space,
    # so we set a dummy objective function (e.g., f(p) = 0).
    def objective_func(p):
        return 0

    # Create the constraints
    # The constraints are of the form 'fun(p) >= 0'
    constraints = []

    # 1. Radial constraints for each sphere
    for i in range(N):
        # Using a closure to capture the value of i at definition time
        def radial_constraint(p, i=i):
            xi, yi = p[3 * i], p[3 * i + 1]
            return (R_c - r_s)**2 - (xi**2 + yi**2)
        constraints.append({'type': 'ineq', 'fun': radial_constraint})

    # 2. Non-overlapping constraints between each pair of spheres
    for i in range(N):
        for j in range(i + 1, N):
            # Using a closure to capture i and j
            def non_overlap_constraint(p, i=i, j=j):
                xi, yi, zi = p[3 * i], p[3 * i + 1], p[3 * i + 2]
                xj, yj, zj = p[3 * j], p[3 * j + 1], p[3 * j + 2]
                return (xi - xj)**2 + (yi - yj)**2 + (zi - zj)**2 - (2 * r_s)**2
            constraints.append({'type': 'ineq', 'fun': non_overlap_constraint})

    # Create bounds for the variables (sphere centers)
    # x is in [-(R_c-r_s), R_c-r_s]
    # y is in [-(R_c-r_s), R_c-r_s]
    # z is in [r_s, H_c-r_s]
    bounds = []
    for i in range(N):
        bounds.extend([
            (-(R_c - r_s), R_c - r_s),
            (-(R_c - r_s), R_c - r_s),
            (r_s, H_c - r_s)
        ])

    # Initial guess: random positions inside the allowed cylinder volume
    p0 = np.random.rand(3 * N)
    for i in range(N):
        p0[3*i] = np.random.uniform(-(R_c - r_s), R_c - r_s) # x
        p0[3*i+1] = np.random.uniform(-(R_c - r_s), R_c - r_s) # y
        p0[3*i+2] = np.random.uniform(r_s, H_c - r_s) # z
        
    print(f"\nAttempting to solve for N = {N} spheres.")
    print("This is a difficult non-convex problem and may take a long time or fail to converge.")

    # Use SLSQP solver. It's not designed for non-convex problems but is a standard choice.
    result = minimize(objective_func, p0, method='SLSQP', bounds=bounds, constraints=constraints, options={'maxiter': 200, 'disp': True})
    
    return result

if __name__ == '__main__':
    # Based on literature review, the most likely answer is 64.
    # We will use this number for the demonstration.
    num_spheres_to_check = 64
    
    result = solve_sphere_packing(num_spheres_to_check)
    
    print("\n--- Solver Results ---")
    print(f"Success: {result.success}")
    print(f"Message: {result.message}")
    if result.success:
        print(f"A feasible arrangement for {num_spheres_to_check} spheres was found.")
    else:
        print(f"Could not find a feasible arrangement for {num_spheres_to_check} spheres (this does not prove it's impossible).")
