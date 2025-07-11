import numpy as np
from scipy.optimize import differential_evolution
import math

# --- Problem Parameters ---
R_CYLINDER = 10.0
H_CYLINDER = 11.0
R_SPHERE = 2.0
D_SPHERE_SQ = (2 * R_SPHERE)**2

# --- Optimization Settings ---
# Note: Higher values for maxiter and popsize give better results but are slower.
MAX_ITERATIONS = 1000
POPULATION_SIZE = 20

def check_feasibility(n_spheres, r_cyl, h_cyl, r_sph):
    """
    Checks if n_spheres can be feasibly packed into the cylinder.
    Returns True if a valid packing is found, False otherwise.
    """
    if n_spheres <= 0:
        return True

    # The penalty function to be minimized.
    # Its value is 0 for a valid packing and > 0 if constraints are violated.
    def penalty(coords):
        coords = coords.reshape(n_spheres, 3)
        total_penalty = 0.0
        
        # 1. Penalty for sphere overlaps
        for i in range(n_spheres):
            for j in range(i + 1, n_spheres):
                dist_sq = np.sum((coords[i] - coords[j])**2)
                if dist_sq < D_SPHERE_SQ:
                    # Penalize based on how much they overlap
                    total_penalty += (D_SPHERE_SQ - dist_sq)

        # 2. Penalty for spheres outside cylinder radius
        max_r_center_sq = (r_cyl - r_sph)**2
        for i in range(n_spheres):
            r_center_sq = coords[i, 0]**2 + coords[i, 1]**2
            if r_center_sq > max_r_center_sq:
                total_penalty += (r_center_sq - max_r_center_sq)
        
        return total_penalty

    # --- Set Bounds for Sphere Center Coordinates ---
    # The height constraints are handled by these bounds.
    max_xy = r_cyl - r_sph
    min_z = r_sph
    max_z = h_cyl - r_sph
    
    bounds = []
    for _ in range(n_spheres):
        bounds.extend([(-max_xy, max_xy), (-max_xy, max_xy), (min_z, max_z)])

    # --- Run the Optimizer ---
    # We use differential_evolution, a global optimizer suitable for such problems.
    result = differential_evolution(
        penalty,
        bounds,
        maxiter=MAX_ITERATIONS, 
        popsize=POPULATION_SIZE,
        tol=1e-5, # The optimization terminates when fun <= tol
        polish=True,
        seed=42 # For reproducible results
    )

    # If the final penalty is very close to zero, we consider the packing feasible.
    is_feasible = result.fun < 1e-6
    print(f"Result for N = {n_spheres}: Feasible = {is_feasible} (Final penalty = {result.fun:.6f})")
    return is_feasible

def find_max_spheres():
    """
    Finds the maximum number of spheres by searching downwards from an upper bound.
    """
    # An upper bound based on volume ratio and a typical packing density (~0.65)
    # This is just a starting point for the search.
    volume_cylinder = math.pi * R_CYLINDER**2 * H_CYLINDER
    volume_sphere = (4/3) * math.pi * R_SPHERE**3
    upper_bound = math.floor((volume_cylinder / volume_sphere) * 0.65)
    
    print(f"Cylinder dimensions: Radius = {R_CYLINDER}, Height = {H_CYLINDER}")
    print(f"Sphere radius: {R_SPHERE}")
    print(f"Starting search for maximum number of spheres from N = {upper_bound} downwards...")
    
    max_n = 0
    for n in range(upper_bound, 0, -1):
        if check_feasibility(n, R_CYLINDER, H_CYLINDER, R_SPHERE):
            max_n = n
            print("\n---------------------------------------------------------")
            print(f"Found the maximum number of spheres: {max_n}")
            print("---------------------------------------------------------")
            # For each number in the final equation (N = 54), print it
            print(f"N = {max_n}")
            break
            
    if max_n == 0:
        print("Could not find a feasible packing for N > 0.")
        
    return max_n

if __name__ == '__main__':
    find_max_spheres()