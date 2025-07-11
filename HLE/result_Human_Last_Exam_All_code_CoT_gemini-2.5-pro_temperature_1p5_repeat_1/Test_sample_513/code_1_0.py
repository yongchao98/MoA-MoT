import numpy as np
from scipy.optimize import minimize

# --- Problem Parameters ---
SPHERE_RADIUS = 2.0
CYLINDER_RADIUS = 10.0
CYLINDER_HEIGHT = 11.0

# --- Derived Parameters ---
# Effective radius for sphere centers
R_EFF = CYLINDER_RADIUS - SPHERE_RADIUS
R_EFF_SQ = R_EFF**2
# Min/max height for sphere centers
Z_MIN = SPHERE_RADIUS
Z_MAX = CYLINDER_HEIGHT - SPHERE_RADIUS
# Minimum squared distance between sphere centers
MIN_DIST_SQ = (2 * SPHERE_RADIUS)**2

# --- Optimizer Settings ---
# Number of random starting configurations to try for each 'n'
NUM_TRIALS = 30
# Max iterations for the optimizer
MAX_ITER = 2000
# Tolerance for considering the solution valid (close to zero violation)
TOLERANCE = 1e-7

def objective_function(x, n):
    """
    Calculates the total violation of packing constraints.
    The goal of the optimizer is to make this function return 0.
    """
    # Reshape the flat input vector 'x' into a list of 3D coordinates
    coords = x.reshape((n, 3))
    
    total_violation = 0.0

    # 1. Non-overlapping constraint violation
    for i in range(n):
        for j in range(i + 1, n):
            dist_sq = np.sum((coords[i] - coords[j])**2)
            # Add penalty if spheres are too close
            violation = max(0, MIN_DIST_SQ - dist_sq)
            total_violation += violation**2 # Use square for smoother gradient

    # 2. Radial constraint violation
    for i in range(n):
        radial_sq = coords[i, 0]**2 + coords[i, 1]**2
        # Add penalty if sphere center is outside the effective radius
        violation = max(0, radial_sq - R_EFF_SQ)
        total_violation += violation**2
        
    return total_violation

def check_packing(n):
    """
    Tries to find a valid packing for 'n' spheres.
    Runs the optimizer multiple times from random starting points.
    """
    print(f"Checking if {n} spheres can be packed...", flush=True)

    # Define the bounds for each coordinate of each sphere center
    # [-R_eff, R_eff] for x and y, [Z_min, Z_max] for z
    bounds = []
    for _ in range(n):
        bounds.extend([(-R_EFF, R_EFF), (-R_EFF, R_EFF), (Z_MIN, Z_MAX)])

    for i in range(NUM_TRIALS):
        # Generate a random initial configuration
        x0 = np.random.rand(n * 3)
        initial_coords = x0.reshape((n, 3))
        initial_coords[:, 0] = initial_coords[:, 0] * 2 * R_EFF - R_EFF
        initial_coords[:, 1] = initial_coords[:, 1] * 2 * R_EFF - R_EFF
        initial_coords[:, 2] = initial_coords[:, 2] * (Z_MAX - Z_MIN) + Z_MIN
        
        # Run the optimizer
        result = minimize(
            objective_function,
            initial_coords.flatten(),
            args=(n,),
            method='L-BFGS-B',
            bounds=bounds,
            options={'maxiter': MAX_ITER, 'ftol': 1e-9}
        )

        # If the final violation is effectively zero, a solution is found
        if result.fun < TOLERANCE:
            print(f"Success! Found a valid packing for {n} spheres.\n")
            return True, result.x.reshape((n, 3))

    print(f"Failed to find a packing for {n} spheres after {NUM_TRIALS} trials.")
    return False, None

def find_max_spheres():
    """
    Searches for the maximum number of spheres that can be packed.
    Starts from a high estimate and decreases n.
    """
    # Start search from a reasonable upper estimate.
    # A loose volume-based estimate suggests n < 80. Let's start at 60.
    for n in range(60, 0, -1):
        is_feasible, final_coords = check_packing(n)
        if is_feasible:
            print("--- OPTIMAL PACKING FOUND ---")
            print(f"Sphere Radius (r): {SPHERE_RADIUS}")
            print(f"Cylinder Radius (R): {CYLINDER_RADIUS}")
            print(f"Cylinder Height (H): {CYLINDER_HEIGHT}")
            print("\nFinal Result:")
            print(f"The maximum number of spheres that can be packed is: {n}")
            
            # Optional: Print the coordinates of the packed spheres
            # print("\nCoordinates of sphere centers (x, y, z):")
            # for i in range(n):
            #    print(f"Sphere {i+1:>2}: ({final_coords[i, 0]:>6.3f}, {final_coords[i, 1]:>6.3f}, {final_coords[i, 2]:>6.3f})")
            
            return n

    return 0

if __name__ == '__main__':
    max_n = find_max_spheres()
    # The required format is to put the answer at the end
    print(f"\n<<<solution>>>")
    print(f"<<<{max_n}>>>")
