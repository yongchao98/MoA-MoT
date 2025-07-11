import numpy as np

def calculate_field_aligned_system(B0):
    """
    Calculates the basis vectors for a field-aligned coordinate system (FAC).
    
    The system is defined as:
    - z_hat is parallel to the local magnetic field B0.
    - y_hat is perpendicular to the plane formed by B0 and a reference vector (here, the radial direction).
    - x_hat completes the right-handed system.

    Args:
        B0 (np.array): The local average magnetic field vector in a base system (e.g., RTN).

    Returns:
        tuple: A tuple containing the three basis vectors (x_hat, y_hat, z_hat).
    """
    # The reference vector is the radial direction, which is the R-axis in an RTN system.
    r_hat = np.array([1.0, 0.0, 0.0])
    
    # The parallel direction (z_hat) is along the local magnetic field B0.
    z_hat = B0 / np.linalg.norm(B0)
    
    # The first perpendicular direction (y_hat) is given by the cross product
    # of the parallel vector and the reference radial vector.
    # This creates a vector perpendicular to the B-R plane.
    y_hat_unnormalized = np.cross(z_hat, r_hat)
    y_hat = y_hat_unnormalized / np.linalg.norm(y_hat_unnormalized)
    
    # The second perpendicular direction (x_hat) completes the right-handed system.
    # x_hat = y_hat x z_hat
    x_hat = np.cross(y_hat, z_hat)
    
    return x_hat, y_hat, z_hat

def transform_to_fac(vector, basis_vectors):
    """
    Transforms a vector into the field-aligned coordinate system.

    Args:
        vector (np.array): The vector to transform (e.g., a magnetic field fluctuation).
        basis_vectors (tuple): The FAC basis vectors (x_hat, y_hat, z_hat).

    Returns:
        np.array: The components of the vector in the FAC.
    """
    x_hat, y_hat, z_hat = basis_vectors
    # Project the vector onto each basis vector to get the new components
    comp_x = np.dot(vector, x_hat)
    comp_y = np.dot(vector, y_hat)
    comp_z = np.dot(vector, z_hat) # This is the parallel component
    return np.array([comp_x, comp_y, comp_z])

# --- Main execution ---

# 1. Define a sample average magnetic field vector B0 at L1 in RTN coordinates (R, T, N).
# This example represents a Parker Spiral field at ~45 degrees.
# B0 = [Br, Bt, Bn] in nanoTeslas (nT)
B0 = np.array([4.0, -4.0, 1.0])

# 2. Define a sample magnetic field fluctuation vector, delta_B, also in RTN.
# This represents the wave itself.
delta_B = np.array([0.1, 0.5, -0.2])

print("--- Inputs in RTN Coordinate System ---")
print(f"Average Magnetic Field B0 (R,T,N): [{B0[0]:.2f}, {B0[1]:.2f}, {B0[2]:.2f}] nT")
print(f"Magnetic Fluctuation delta_B (R,T,N): [{delta_B[0]:.2f}, {delta_B[1]:.2f}, {delta_B[2]:.2f}] nT\n")

# 3. Calculate the basis vectors for the Field-Aligned Coordinate (FAC) system.
# The components of the basis vectors are expressed in the original RTN system.
fac_x_hat, fac_y_hat, fac_z_hat = calculate_field_aligned_system(B0)

print("--- Field-Aligned Coordinate (FAC) Basis Vectors (in RTN) ---")
print(f"Parallel direction (z_hat):      [{fac_z_hat[0]:.3f}, {fac_z_hat[1]:.3f}, {fac_z_hat[2]:.3f}]")
print(f"Perpendicular direction 1 (x_hat): [{fac_x_hat[0]:.3f}, {fac_x_hat[1]:.3f}, {fac_x_hat[2]:.3f}]")
print(f"Perpendicular direction 2 (y_hat): [{fac_y_hat[0]:.3f}, {fac_y_hat[1]:.3f}, {fac_y_hat[2]:.3f}]\n")


# 4. Transform the fluctuation vector delta_B into the new FAC system.
delta_B_fac = transform_to_fac(delta_B, (fac_x_hat, fac_y_hat, fac_z_hat))

print("--- Final Transformed Fluctuation Vector ---")
print("The helicity calculation should use the two perpendicular components (B_perp1, B_perp2).")
print(f"Component Perpendicular 1 (along x_hat): {delta_B_fac[0]:.4f} nT")
print(f"Component Perpendicular 2 (along y_hat): {delta_B_fac[1]:.4f} nT")
print(f"Component Parallel (along z_hat):        {delta_B_fac[2]:.4f} nT")

# We can express the final result as an equation:
# delta_B_fac = [B_perp1, B_perp2, B_parallel]
print("\nFinal Equation: delta_B_fac = [{:.4f}, {:.4f}, {:.4f}]".format(delta_B_fac[0], delta_B_fac[1], delta_B_fac[2]))
