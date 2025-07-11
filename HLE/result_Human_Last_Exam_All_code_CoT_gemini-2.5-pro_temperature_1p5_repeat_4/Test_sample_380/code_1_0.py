import numpy as np

def calculate_field_aligned_components():
    """
    This function demonstrates the correct procedure for finding the magnetic field
    components perpendicular to the local mean magnetic field.
    """
    # 1. Define the Mean Magnetic Field (B0) in a standard frame (e.g., GSE).
    # This vector represents a typical Parker Spiral field at L1, not purely radial.
    # The X component is radial, Y is tangential. Units are nanoTeslas (nT).
    # A purely radial field would be [-5, 0, 0].
    B0 = np.array([-4.0, 4.0, 1.0])

    # 2. Define a sample instantaneous magnetic field fluctuation (delta_B).
    # This is the wave's magnetic field vector that we want to analyze.
    delta_B = np.array([0.5, 0.8, -0.2])

    print("--- Input Vectors (in GSE frame) ---")
    print(f"Mean Magnetic Field B0 = {B0} nT")
    print(f"Fluctuating Field delta_B = {delta_B} nT")
    print("-" * 35)

    # 3. Create the Field-Aligned Coordinate (FAC) system based on B0.
    
    # The parallel unit vector (z') is in the direction of B0.
    B0_magnitude = np.linalg.norm(B0)
    z_prime_axis = B0 / B0_magnitude

    # To define the perpendicular plane, we need another reference vector.
    # We can use the GSE Z-axis, as long as it's not parallel to B0.
    ref_vec = np.array([0.0, 0.0, 1.0])
    
    # The y' axis is perpendicular to both B0 and the reference vector.
    y_prime_axis = np.cross(z_prime_axis, ref_vec)
    y_prime_axis = y_prime_axis / np.linalg.norm(y_prime_axis)

    # The x' axis completes the right-handed system (x' = y' cross z').
    x_prime_axis = np.cross(y_prime_axis, z_prime_axis)

    print("\n--- Field-Aligned Coordinate (FAC) System Basis Vectors ---")
    # Using f-string formatting to show each number clearly
    print(f"Parallel axis (z'):    [{z_prime_axis[0]:.4f}, {z_prime_axis[1]:.4f}, {z_prime_axis[2]:.4f}]")
    print(f"Perp. axis 1 (x'):    [{x_prime_axis[0]:.4f}, {x_prime_axis[1]:.4f}, {x_prime_axis[2]:.4f}]")
    print(f"Perp. axis 2 (y'):    [{y_prime_axis[0]:.4f}, {y_prime_axis[1]:.4f}, {y_prime_axis[2]:.4f}]")
    print("-" * 35)

    # 4. Project the fluctuation delta_B onto the new FAC axes.
    # These dot products give the components of delta_B in the new system.
    
    delta_B_parallel = np.dot(delta_B, z_prime_axis)
    delta_B_perp_x = np.dot(delta_B, x_prime_axis)
    delta_B_perp_y = np.dot(delta_B, y_prime_axis)

    # The helicity calculation would use delta_B_perp_x and delta_B_perp_y.
    
    print("\n--- Fluctuation Components in FAC System ---")
    print("These are the physically meaningful components for wave analysis.")
    print(f"Component Parallel to B0:       {delta_B_parallel:.4f} nT")
    print(f"Component Perpendicular to B0 (x'): {delta_B_perp_x:.4f} nT")
    print(f"Component Perpendicular to B0 (y'): {delta_B_perp_y:.4f} nT")
    print("-" * 35)

    print("\nConclusion: To calculate magnetic helicity for an AIC wave, you must use")
    print("the components perpendicular to the local magnetic field (e.g., the x' and y'")
    print(f"components calculated above), which are {delta_B_perp_x:.4f} and {delta_B_perp_y:.4f} nT in this example.")


if __name__ == '__main__':
    calculate_field_aligned_components()
