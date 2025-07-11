import numpy as np
from scipy.optimize import minimize

# --- Problem Parameters ---
R_CYLINDER = 10.0  # Cylinder radius
H_CYLINDER = 11.0  # Cylinder height
R_SPHERE = 2.0     # Sphere radius

# --- Derived Constraints for Sphere Centers ---
# Effective radius for centers: cylinder radius minus sphere radius
R_EFFECTIVE = R_CYLINDER - R_SPHERE
R_EFFECTIVE_SQ = R_EFFECTIVE**2
# Min and max height for centers
Z_MIN = R_SPHERE
Z_MAX = H_CYLINDER - R_SPHERE
# Squared diameter of spheres for non-overlap constraint
DIAMETER_SQ = (2 * R_SPHERE)**2


def objective_function(positions, n_spheres):
    """
    Calculates the penalty score for a given arrangement of sphere centers.
    A score of 0 means a valid packing.
    """
    # Reshape the 1D array of positions into a (n_spheres, 3) array
    positions = positions.reshape((n_spheres, 3))
    
    total_penalty = 0.0

    # 1. Penalty for spheres overlapping each other
    # Calculate pairwise squared distances in a vectorized way
    p_tiled = np.tile(positions, (n_spheres, 1, 1))
    p_transposed = np.transpose(p_tiled, axes=(1, 0, 2))
    dist_sq_matrix = np.sum((p_tiled - p_transposed)**2, axis=2)
    
    # We only need to check the upper triangle of the distance matrix
    # for unique pairs (i, j) where i < j.
    indices = np.triu_indices(n_spheres, k=1)
    dist_sq_pairs = dist_sq_matrix[indices]
    
    # Calculate penalty for any pair that is too close
    overlap_dist = DIAMETER_SQ - dist_sq_pairs
    overlap_penalty = np.sum(np.maximum(0, overlap_dist)**2)
    total_penalty += overlap_penalty

    # 2. Penalty for spheres outside the cylinder radius
    # The Z-bounds are handled by the optimizer's 'bounds' parameter.
    radial_pos_sq = np.sum(positions[:, :2]**2, axis=1)
    outside_dist = radial_pos_sq - R_EFFECTIVE_SQ
    radial_penalty = np.sum(np.maximum(0, outside_dist)**2)
    total_penalty += radial_penalty
    
    return total_penalty

def find_packing(n_spheres, num_trials=20):
    """
    Tries to find a valid packing for n_spheres.
    Returns the positions if successful, otherwise None.
    """
    print(f"Trying to pack N = {n_spheres} spheres...")

    # Define the bounds for each variable (x, y, z for each sphere)
    bounds = []
    for _ in range(n_spheres):
        bounds.extend([
            (-R_EFFECTIVE, R_EFFECTIVE),
            (-R_EFFECTIVE, R_EFFECTIVE),
            (Z_MIN, Z_MAX)
        ])

    for i in range(num_trials):
        # Generate a random initial guess for the sphere positions
        initial_guess = np.zeros(3 * n_spheres)
        for k in range(n_spheres):
            # Use rejection sampling to get a random point within the cylinder's base circle
            while True:
                x = np.random.uniform(-R_EFFECTIVE, R_EFFECTIVE)
                y = np.random.uniform(-R_EFFECTIVE, R_EFFECTIVE)
                if x**2 + y**2 <= R_EFFECTIVE_SQ:
                    break
            z = np.random.uniform(Z_MIN, Z_MAX)
            initial_guess[3*k : 3*k+3] = [x, y, z]

        # Run the optimizer
        res = minimize(
            fun=objective_function,
            x0=initial_guess,
            args=(n_spheres,),
            method='L-BFGS-B',
            bounds=bounds,
            options={'ftol': 1e-9, 'maxiter': 500}
        )

        # If the final penalty is very close to zero, we found a valid packing
        if res.success and res.fun < 1e-6:
            print(f"Success! Found a valid packing for N = {n_spheres} on trial {i+1}.")
            return res.x.reshape((n_spheres, 3))

    print("Failed to find a packing after all trials.")
    return None

def solve_packing_problem():
    """
    Searches for the maximum number of spheres that can be packed.
    """
    # Search downwards from a reasonable upper estimate
    for n in range(70, 50, -1):
        positions = find_packing(n)
        if positions is not None:
            print("\n" + "="*50)
            print("                OPTIMAL PACKING FOUND")
            print("="*50)
            print(f"Sphere radius (r): {R_SPHERE}")
            print(f"Cylinder radius (R): {R_CYLINDER}")
            print(f"Cylinder height (H): {H_CYLINDER}")
            print("\nThe final equation for the maximum number of spheres is:")
            print(f"N = {n}")
            print("="*50)
            return n
    
    print("Could not find a solution in the searched range.")
    return 0

if __name__ == '__main__':
    max_spheres = solve_packing_problem()
    # The final answer is submitted below based on the result of this method.
    # Note: The script is non-deterministic and runtime may vary. 
    # Increasing `num_trials` in find_packing increases the chance of success.