import numpy as np
from scipy.optimize import differential_evolution
import time

# --- Problem Parameters ---
R_CYLINDER = 10.0  # Cylinder radius: 10cm
H_CYLINDER = 11.0  # Cylinder height: 11cm
R_SPHERE = 2.0     # Sphere radius: 2cm

# --- Derived Parameters for Constraints ---
# The effective squared radius for the sphere centers
R_EFFECTIVE_SQ = (R_CYLINDER - R_SPHERE)**2  
# Min and max height for sphere centers
H_MIN_CENTER = R_SPHERE
H_MAX_CENTER = H_CYLINDER - R_SPHERE
# The squared minimum distance between any two sphere centers (diameter squared)
MIN_DIST_SQ = (2 * R_SPHERE)**2

def objective_function(positions, n_spheres):
    """
    Calculates a penalty score for a given arrangement of spheres.
    The goal of the optimizer is to find an arrangement where this score is zero.
    A score of zero means all constraints are satisfied (a valid packing).
    
    Args:
        positions (np.array): A flat array of coordinates [x1, y1, z1, x2, y2, z2, ...].
        n_spheres (int): The number of spheres in the current test.

    Returns:
        float: The total penalty score. A score of 0 indicates success.
    """
    # Reshape the flat array into a (n_spheres, 3) array of center coordinates
    centers = positions.reshape((n_spheres, 3))
    
    total_penalty = 0.0
    
    # 1. Containment Penalty: For spheres outside the cylinder.
    # a) Radial constraint: x^2 + y^2 <= (R-r)^2
    radial_dist_sq = np.sum(centers[:, :2]**2, axis=1)
    total_penalty += np.sum(np.maximum(0, radial_dist_sq - R_EFFECTIVE_SQ))
    
    # b) Height constraint: r <= z <= H-r
    z_coords = centers[:, 2]
    total_penalty += np.sum(np.maximum(0, H_MIN_CENTER - z_coords))
    total_penalty += np.sum(np.maximum(0, z_coords - H_MAX_CENTER))

    # 2. Overlap Penalty: For any two spheres that are too close.
    if n_spheres > 1:
        # Vectorized calculation of all pairwise squared distances
        pdist_sq = np.sum((centers[:, np.newaxis, :] - centers[np.newaxis, :, :])**2, axis=-1)
        
        # We only need the upper triangle of the distance matrix to count each pair once
        indices = np.triu_indices(n_spheres, k=1)
        pairwise_sq_dists = pdist_sq[indices]
        
        # Add penalty for pairs where distance squared is less than minimum distance squared
        total_penalty += np.sum(np.maximum(0, MIN_DIST_SQ - pairwise_sq_dists))
        
    return total_penalty

def find_max_spheres():
    """
    Iteratively searches for the maximum number of spheres that can be packed.
    It starts from a high estimate and decrements until a valid packing is found.
    """
    # Based on heuristic analysis, the solution is likely in the 50-60 range.
    # We will search downwards from 60.
    for n in range(60, 50, -1):
        print(f"\n--- Checking for N = {n} spheres ---")
        start_time = time.time()
        
        # Define the search space (bounds) for each variable (x_i, y_i, z_i).
        # The circular constraint is handled by the penalty function, so we use a square box here.
        bounds = []
        for _ in range(n):
            bounds.extend([
                (-(R_CYLINDER - R_SPHERE), R_CYLINDER - R_SPHERE),  # x_i
                (-(R_CYLINDER - R_SPHERE), R_CYLINDER - R_SPHERE),  # y_i
                (H_MIN_CENTER, H_MAX_CENTER)                      # z_i
            ])
            
        # Use SciPy's differential_evolution, a global optimization algorithm.
        result = differential_evolution(
            func=objective_function,
            bounds=bounds,
            args=(n,),
            strategy='best1bin',
            maxiter=1000,   # Iterations for a thorough search
            popsize=30,     # Larger population helps find global minimum
            tol=1e-8,       # Tolerance for convergence
            seed=0,         # Use a seed for reproducible results
            workers=-1      # Use all available CPU cores to speed up
        )
        
        end_time = time.time()
        print(f"Time taken for N = {n}: {end_time - start_time:.2f} seconds")
        print(f"Final penalty score for N = {n}: {result.fun:.6f}")

        # If the penalty is effectively zero, we have found a valid packing.
        if result.fun < 1e-6:
            print("\n=======================================================")
            print(f"Found a feasible packing for N = {n} spheres.")
            print(f"This is the maximum number of spheres that can be optimally packed.")
            print(f"Cylinder (Radius={R_CYLINDER}, Height={H_CYLINDER}), Sphere (Radius={R_SPHERE})")
            print(f"The maximum number of spheres is: {n}")
            print("=======================================================")
            return n
            
    print("\nCould not find a solution in the specified range with the given parameters.")
    return None

if __name__ == '__main__':
    max_spheres = find_max_spheres()
    if max_spheres is None:
        print("Optimization failed to find a solution. Try increasing 'maxiter' or 'popsize'.")
