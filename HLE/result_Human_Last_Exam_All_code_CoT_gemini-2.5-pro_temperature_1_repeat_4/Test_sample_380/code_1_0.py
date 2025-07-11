import numpy as np

def transform_to_field_aligned_coords(mean_field, fluctuation_vector):
    """
    Transforms a fluctuation vector to a Field-Aligned Coordinate (FAC) system.

    In this FAC system:
    - z_fac is parallel to the mean magnetic field.
    - y_fac is perpendicular to the plane formed by the mean field and the original Z-axis.
    - x_fac completes the right-handed system.

    Args:
        mean_field (np.ndarray): A 3-element array representing the mean magnetic field vector [Bx, By, Bz].
        fluctuation_vector (np.ndarray): A 3-element array for the fluctuation vector [dBx, dBy, dBz].

    Returns:
        np.ndarray: The fluctuation vector transformed into the FAC system.
    """
    # 1. Define the basis vectors of the new FAC system
    
    # The z'-axis is parallel to the mean magnetic field B0
    z_fac = mean_field / np.linalg.norm(mean_field)
    
    # The y'-axis is perpendicular to the plane containing B0 and the original Z-axis ([0, 0, 1])
    # This is a common convention to define the perpendicular plane.
    original_z_axis = np.array([0., 0., 1.])
    y_fac = np.cross(z_fac, original_z_axis)
    y_fac = y_fac / np.linalg.norm(y_fac)
    
    # The x'-axis completes the right-handed system (x' = y' x z')
    x_fac = np.cross(y_fac, z_fac)

    # 2. Create the rotation matrix to transform from the original system to the FAC system.
    # The rows of the rotation matrix are the new basis vectors.
    rotation_matrix = np.array([x_fac, y_fac, z_fac])
    
    # 3. Apply the rotation to the fluctuation vector
    fluctuation_in_fac = rotation_matrix.dot(fluctuation_vector)
    
    return fluctuation_in_fac, rotation_matrix

# --- Example Usage ---

# Assume our original data is in a radial-based system (e.g., RTN)
# R is [1, 0, 0], T is [0, 1, 0], N is [0, 0, 1]

# Define a mean magnetic field (B0) typical of the Parker Spiral at 1 AU.
# It has both radial (Br) and tangential (Bt) components. Not purely radial.
# Let's say B0 = [4 nT, 4 nT, 1 nT] in the RTN system.
B0 = np.array([4.0, 4.0, 1.0])

# Define a sample magnetic field fluctuation vector (delta_B)
delta_B_original = np.array([0.5, -0.2, 0.8])

# Perform the transformation
delta_B_fac, rot_matrix = transform_to_field_aligned_coords(B0, delta_B_original)

# The components used for helicity calculation are the first two (perpendicular)
# components of the transformed vector.
delta_B_perp_1 = delta_B_fac[0]
delta_B_perp_2 = delta_B_fac[1]
delta_B_parallel = delta_B_fac[2]

# --- Print the results ---
print("--- Justification for using a Field-Aligned Coordinate System ---")
print(f"Original Mean Magnetic Field (B0) in RTN system: {B0}")
print("Note: This B0 is not purely radial, confirming the Parker Spiral model.\n")

print(f"Original Fluctuation Vector (delta_B) in RTN system: {delta_B_original}\n")

print("Rotation Matrix (from RTN to FAC):")
print(rot_matrix)
print("\nThis matrix rotates vectors so the new z-axis is parallel to B0.\n")

print(f"Fluctuation Vector in Field-Aligned Coordinates (delta_B_fac): {delta_B_fac}\n")

print("--- Components for Helicity Calculation ---")
print(f"Component 1 perpendicular to B0: {delta_B_perp_1:.4f}")
print(f"Component 2 perpendicular to B0: {delta_B_perp_2:.4f}")
print(f"Component parallel to B0: {delta_B_parallel:.4f}\n")
print("The helicity calculation for AIC waves correctly uses the two perpendicular components shown above.")
