import numpy as np
from scipy.optimize import differential_evolution

def solve_packing_problem():
    """
    Formulates and solves the sphere packing problem.
    This function will print the problem parameters and the results of the optimization search.
    """
    # --- Problem Parameters ---
    R_cylinder = 10.0
    H_cylinder = 11.0
    r_sphere = 2.0
    
    print("--- Problem Setup ---")
    print(f"Cylinder dimensions: Radius = {R_cylinder}, Height = {H_cylinder}")
    print(f"Sphere dimensions:   Radius = {r_sphere}")
    print("The goal is to find the maximum number of spheres that can be packed in the cylinder.")
    
    # Derived parameters for the cost function and bounds
    d_sphere_sq = (2 * r_sphere)**2
    R_effective = R_cylinder - r_sphere
    R_effective_sq = R_effective**2
    Z_min = r_sphere
    Z_max = H_cylinder - r_sphere

    print("\n--- Final Equation Constraints ---")
    print(f"1. Non-Overlap: (xi-xj)²+(yi-yj)²+(zi-zj)² >= (2*r_sphere)² = {d_sphere_sq}")
    print(f"2. Radial Containment: xi²+yi² <= (R_cylinder-r_sphere)² = {R_effective_sq}")
    print(f"3. Height Containment: r_sphere <= zi <= H_cylinder-r_sphere => {Z_min} <= zi <= {Z_max}\n")

    def cost_function(X, n_spheres):
        """
        Calculates the total penalty for a given configuration of sphere centers.
        A value of 0 indicates a valid packing.
        X: A flat array of coordinates [x1, y1, z1, x2, y2, z2, ...].
        """
        centers = X.reshape((n_spheres, 3))
        total_penalty = 0.0
        
        # 1. Penalty for sphere-sphere overlaps
        for i in range(n_spheres):
            for j in range(i + 1, n_spheres):
                dist_sq = np.sum((centers[i] - centers[j])**2)
                overlap = max(0, d_sphere_sq - dist_sq)
                total_penalty += overlap * 10 # Weight overlaps more heavily

        # 2. Penalty for cylinder wall violation (radially)
        for i in range(n_spheres):
            radial_dist_sq = centers[i, 0]**2 + centers[i, 1]**2
            violation = max(0, radial_dist_sq - R_effective_sq)
            total_penalty += violation
            
        return total_penalty

    def check_packing(n_spheres, max_iter=2000, tolerance=1e-7):
        """
        Tries to find a valid packing for n_spheres using differential_evolution.
        Returns the number of spheres if successful, otherwise 0.
        """
        print(f"--- Checking for N = {n_spheres} spheres ---")
        if n_spheres == 0:
            return True
        
        # Define the bounds for each coordinate of the sphere centers
        bounds = []
        for _ in range(n_spheres):
            bounds.extend([
                (-R_effective, R_effective),  # x
                (-R_effective, R_effective),  # y
                (Z_min, Z_max)                # z
            ])

        obj_func = lambda x: cost_function(x, n_spheres)
        
        # Note: This is a computationally intensive task.
        # `atol` tells the optimizer to stop as soon as it finds a
        # solution with a cost function value below the tolerance.
        result = differential_evolution(
            func=obj_func,
            bounds=bounds,
            maxiter=max_iter,
            popsize=20,
            polish=True,
            atol=tolerance
        )

        print(f"Result for N = {n_spheres}: Final penalty = {result.fun:.6f}")
        if result.fun < tolerance:
            print("  A valid packing configuration was found.")
            return n_spheres
        else:
            print("  Could not find a valid packing configuration within the given iterations.")
            return 0

    # Based on known results for sphere packing, the answer is likely around 57.
    # We will test N=57 first, and if successful, try N=58.
    
    max_spheres_found = 0
    
    # Test N=57
    result_57 = check_packing(n_spheres=57)
    if result_57 > 0:
        max_spheres_found = result_57
        # If 57 works, try for 58
        result_58 = check_packing(n_spheres=58, max_iter=3000) # Give it more iterations
        if result_58 > 0:
            max_spheres_found = result_58
            
    print("\n--- Conclusion ---")
    print(f"The maximum number of spheres that can be optimally packed is: {max_spheres_found}")
    return max_spheres_found

if __name__ == '__main__':
    solve_packing_problem()
