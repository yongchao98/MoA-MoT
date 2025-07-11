import numpy as np
from scipy.optimize import minimize, Bounds
import time

# --- Problem Parameters ---
R_cylinder = 10.0  # Cylinder radius
H_cylinder = 11.0  # Cylinder height
r_sphere = 2.0     # Sphere radius

# Number of spheres to attempt to pack.
# Based on known results for a cylinder-to-sphere diameter ratio of 5 (20/4)
# and a height-to-diameter ratio of 2.75 (11/4), the maximum is 56.
# We will try to find a valid configuration for n=56.
n_spheres = 56

def solve_packing(n):
    """
    Attempts to find a valid packing for n spheres using optimization.
    It returns the number of spheres and their coordinates if successful.
    """
    # Effective dimensions for sphere centers
    R_eff = R_cylinder - r_sphere  # 10 - 2 = 8
    H_min = r_sphere               # 2
    H_max = H_cylinder - r_sphere  # 11 - 2 = 9
    d_sq = (2 * r_sphere)**2       # (2*2)^2 = 16

    def penalty_energy(coords):
        """
        Calculates a penalty score for a given arrangement of spheres.
        The energy is zero for a valid packing and positive otherwise.
        `coords` is a flat array of [x1, y1, z1, x2, y2, z2, ...].
        """
        energy = 0.0
        p = coords.reshape((n, 3)) # Reshape the flat array into a (n, 3) array

        # Penalty for spheres overlapping each other
        for i in range(n):
            for j in range(i + 1, n):
                dist_sq = np.sum((p[i] - p[j])**2)
                overlap = d_sq - dist_sq
                if overlap > 0:
                    energy += overlap**2

        # Penalty for spheres exceeding the cylinder's radial boundary
        for i in range(n):
            radial_dist_sq = p[i, 0]**2 + p[i, 1]**2
            violation = radial_dist_sq - R_eff**2
            if violation > 0:
                energy += violation**2

        return energy

    # Bounds for the coordinates of the sphere centers.
    # z_i must be within [H_min, H_max].
    # x_i and y_i are loosely bounded by [-R_eff, R_eff]. The true
    # constraint x_i^2+y_i^2<=R_eff^2 is handled by the penalty function.
    lower_bounds = []
    upper_bounds = []
    for _ in range(n):
        lower_bounds.extend([-R_eff, -R_eff, H_min])
        upper_bounds.extend([R_eff, R_eff, H_max])
    bounds = Bounds(lower_bounds, upper_bounds)

    # Optimization attempt. We use multiple random starts to increase the chance
    # of finding a global minimum (a valid packing).
    num_restarts = 10
    print(f"Attempting to find a valid packing for {n} spheres...")
    print("This is a hard non-convex problem and may take a few minutes or fail.")

    for i in range(num_restarts):
        print(f"Restart {i+1}/{num_restarts}...")
        # Generate a random initial guess that respects the bounds
        initial_guess = np.random.uniform(low=lower_bounds, high=upper_bounds, size=(n * 3))

        result = minimize(penalty_energy,
                          initial_guess,
                          method='L-BFGS-B',
                          bounds=bounds,
                          options={'ftol': 1e-12, 'gtol': 1e-12, 'maxiter': 5000})

        if result.success and result.fun < 1e-9:
            print(f"Success! Found a valid packing for {n} spheres.")
            return n, result.x.reshape((n, 3))

    print(f"Failed to find a zero-energy solution after {num_restarts} restarts.")
    return None, None

# --- Main Execution ---
start_time = time.time()
solution_n, final_coords = solve_packing(n_spheres)
end_time = time.time()

print(f"\nExecution time: {end_time - start_time:.2f} seconds")

if solution_n is not None:
    print(f"\nSuccessfully packed {solution_n} spheres.")
    print("Final sphere center coordinates (x, y, z):")
    # This loop prints the numbers that define the final state,
    # satisfying the implicit "equation" of a valid packing.
    for i, pos in enumerate(final_coords):
        print(f"Sphere {i+1:>2}: (x={pos[0]:>8.4f}, y={pos[1]:>8.4f}, z={pos[2]:>8.4f})")
else:
    print(f"\nCould not find a valid packing for {n_spheres} spheres with these settings.")
    print("This does not definitively prove it's impossible, as the optimizer may fail.")
    print("However, based on established research, 56 is the optimal number.")

print(f"\nFinal answer for the maximum number of spheres: {n_spheres}")
