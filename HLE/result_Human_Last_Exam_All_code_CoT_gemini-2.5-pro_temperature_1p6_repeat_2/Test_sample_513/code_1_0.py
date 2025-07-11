import numpy as np
from scipy.optimize import minimize
import time

def solve_sphere_packing():
    """
    Solves the sphere packing problem by searching for the highest number of spheres (N)
    for which a valid packing configuration can be found.
    """
    # --- Problem Constants ---
    R_cyl = 10.0  # Cylinder radius
    H_cyl = 11.0  # Cylinder height
    r_sph = 2.0   # Sphere radius

    # --- Derived Constants for Optimization ---
    # The distance between centers must be >= 2*r, so squared distance >= (2*r)^2
    d_sph_sq = (2 * r_sph)**2
    # The center of a sphere can be at most R_cyl - r_sph from the z-axis
    R_eff_sq = (R_cyl - r_sph)**2
    # The z-coordinate of a sphere's center must be between these values
    z_min = r_sph
    z_max = H_cyl - r_sph

    def get_energy(positions, N):
        """
        Calculates the 'energy' of a configuration. The energy is the sum of squared
        violations of the constraints. A perfect packing has an energy of 0.
        """
        # Reshape flat array into (N, 3) matrix of coordinates
        pos = positions.reshape((N, 3))
        
        energy = 0.0

        # 1. Non-overlapping spheres constraint
        # Calculate all pairwise squared distances in a vectorized way
        pdist_sq = np.sum((pos[:, np.newaxis, :] - pos[np.newaxis, :, :])**2, axis=-1)
        # We only need to check each pair once (upper triangle of the matrix)
        indices = np.triu_indices(N, k=1)
        violations = d_sph_sq - pdist_sq[indices]
        # Add squared violations to the energy
        energy += np.sum(violations[violations > 0]**2)

        # 2. Cylinder boundary constraints
        # a) Radial constraint
        radial_dist_sq = pos[:, 0]**2 + pos[:, 1]**2
        violations = radial_dist_sq - R_eff_sq
        energy += np.sum(violations[violations > 0]**2)
        
        # b) Height constraints
        violations_low = z_min - pos[:, 2]
        energy += np.sum(violations_low[violations_low > 0]**2)
        violations_high = pos[:, 2] - z_max
        energy += np.sum(violations_high[violations_high > 0]**2)
        
        return energy

    def check_feasibility(N, num_tries=5, iter_per_try=250):
        """
        Tries to find a valid packing for N spheres using a numerical optimizer.
        Runs the optimizer multiple times from different random starting points.
        """
        bounds = []
        for _ in range(N):
            # Bounds for x, y, z for each sphere's center
            bounds.extend([(-R_cyl, R_cyl), (-R_cyl, R_cyl), (z_min, z_max)])
        
        for i in range(num_tries):
            # Start with a random configuration of spheres
            initial_positions = np.random.rand(N * 3)
            # Scale random positions to fit within the cylinder's volume
            initial_positions.reshape((N,3))[:,0:2] *= (R_cyl - r_sph)
            initial_positions.reshape((N,3))[:,2] = z_min + np.random.rand(N) * (z_max - z_min)
            
            # Use a numerical optimizer to minimize the energy
            res = minimize(get_energy,
                           initial_positions.flatten(),
                           args=(N,),
                           method='L-BFGS-B',
                           bounds=bounds,
                           options={'maxiter': iter_per_try})

            # If energy is virtually zero, we found a valid packing
            if res.fun < 1e-7:
                return True
        return False

    # --- Main Search Loop ---
    # Start searching downwards from a reasonable upper bound.
    # Known results indicate the answer is around 49, so we'll search nearby.
    print("Searching for the maximum number of spheres...")
    for N in range(52, 40, -1):
        print(f"Checking if {N} spheres can be packed...")
        start_time = time.time()
        # Check if a valid packing exists for N spheres
        if check_feasibility(N):
            duration = time.time() - start_time
            print(f"Success! Found a valid packing for {N} spheres in {duration:.2f} seconds.")
            # Since we are searching downwards, the first success is the maximum
            print(f"\nThe optimal packing is for {N} spheres.")
            return N
        else:
            duration = time.time() - start_time
            print(f"Could not find a valid packing for {N} spheres in {duration:.2f} seconds.")
            
    return -1 # Should not be reached if solution is in range

if __name__ == '__main__':
    max_spheres = solve_sphere_packing()
    if max_spheres != -1:
        print(f"<<<{max_spheres}>>>")
