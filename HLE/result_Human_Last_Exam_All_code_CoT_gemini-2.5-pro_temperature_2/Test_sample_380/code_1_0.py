import numpy as np

def calculate_perpendicular_components(B_avg, b_fluc):
    """
    Calculates the magnetic field fluctuation components perpendicular to the local
    average magnetic field direction.

    Args:
        B_avg (np.array): A 3-element array representing the average background
                          magnetic field vector (e.g., [Bx, By, Bz]).
        b_fluc (np.array): A 3-element array representing the fluctuating wave
                           magnetic field vector at a moment in time.
    """
    print(f"Average Local Magnetic Field (B_avg): {B_avg}")
    print(f"Fluctuating Magnetic Field (b_fluc): {b_fluc}\n")

    # Step 1: Create a Field-Aligned Coordinate (FAC) system.
    # The first basis vector is parallel to the average magnetic field.
    e_parallel = B_avg / np.linalg.norm(B_avg)

    # To create the perpendicular vectors, we use a fixed reference (here, Z-axis)
    # and the cross product to ensure orthogonality. This defines a right-handed system.
    # Note: If B_avg is perfectly aligned with the reference, another reference must be chosen.
    e_perp2_unnormalized = np.cross(e_parallel, [0, 0, 1])
    e_perp2 = e_perp2_unnormalized / np.linalg.norm(e_perp2_unnormalized)

    e_perp1 = np.cross(e_perp2, e_parallel)

    print("--- Field-Aligned Coordinate System Basis Vectors ---")
    print(f"e_parallel (along B_avg): {np.round(e_parallel, 3)}")
    print(f"e_perp1               : {np.round(e_perp1, 3)}")
    print(f"e_perp2               : {np.round(e_perp2, 3)}\n")

    # Step 2: Project the fluctuating field onto the new basis vectors.
    b_parallel_comp = np.dot(b_fluc, e_parallel)
    b_perp1_comp = np.dot(b_fluc, e_perp1)
    b_perp2_comp = np.dot(b_fluc, e_perp2)

    print("--- Fluctuating Field in Field-Aligned Coordinates ---")
    print(f"Component parallel to B_avg: {b_parallel_comp:.3f}")
    print(f"First component perpendicular to B_avg (b_perp1): {b_perp1_comp:.3f}")
    print(f"Second component perpendicular to B_avg (b_perp2): {b_perp2_comp:.3f}\n")

    # Step 3: The normalized magnetic helicity calculation would then use the time series
    # of b_perp1 and b_perp2. The final equation for helicity involves these components.
    print("Final components needed for the helicity equation are:")
    print(f"b_perp1 = {b_perp1_comp:.3f}")
    print(f"b_perp2 = {b_perp2_comp:.3f}")


# --- Example Usage ---
# At 1 AU, the Parker spiral angle is ~45 degrees. So the field is not radial.
# Let's represent the radial direction as the X-axis. A non-radial field
# would have significant X and Y components.
# B_avg units are in nT.
B_avg_parker_spiral = np.array([-4.0, 4.0, 1.0]) # Example non-radial field

# A fluctuation from a transverse wave.
b_fluctuation = np.array([0.5, 0.7, -1.2])

calculate_perpendicular_components(B_avg_parker_spiral, b_fluctuation)
