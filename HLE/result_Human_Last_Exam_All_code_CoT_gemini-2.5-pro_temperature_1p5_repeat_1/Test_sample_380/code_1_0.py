import numpy as np

def get_perpendicular_fluctuation(B_background, dB_fluctuation):
    """
    Calculates the component of a magnetic field fluctuation that is
    perpendicular to the background magnetic field.

    Args:
        B_background (np.ndarray): A 3-element array representing the background
                                   magnetic field vector [Bx, By, Bz].
        dB_fluctuation (np.ndarray): A 3-element array representing the fluctuation
                                     vector [dBx, dBy, dBz].

    Returns:
        np.ndarray: A 3-element array representing the perpendicular component
                    of the fluctuation vector. The components of this vector
                    are the ones to be used for helicity calculation.
    """
    # Ensure inputs are numpy arrays for vector operations
    B_background = np.array(B_background)
    dB_fluctuation = np.array(dB_fluctuation)

    # 1. Get the unit vector in the direction of the background field (B_parallel_hat)
    norm_B_background = np.linalg.norm(B_background)
    if norm_B_background == 0:
        raise ValueError("Background magnetic field cannot be a zero vector.")
    b_parallel_unit = B_background / norm_B_background

    # 2. Project the fluctuation vector onto the parallel unit vector to find the
    #    magnitude of the parallel component of the fluctuation.
    dB_parallel_magnitude = np.dot(dB_fluctuation, b_parallel_unit)

    # 3. Construct the vector of the parallel component of the fluctuation.
    dB_parallel_vector = dB_parallel_magnitude * b_parallel_unit

    # 4. The perpendicular component is the total fluctuation minus the parallel component.
    dB_perp_vector = dB_fluctuation - dB_parallel_vector

    return dB_perp_vector

# --- Example Usage ---
# Let's assume a Parker Spiral field at L1, mostly in the X-Y plane.
# B0x is radial, B0y is azimuthal. B0z is small (out of ecliptic).
# Units are arbitrary (e.g., nanoTesla, nT).
B0 = [3.5, -3.5, 0.5]

# Now, let's assume we measure a magnetic field fluctuation (the wave).
dB = [0.1, 0.8, -0.4]

# Calculate the perpendicular component of the fluctuation
dB_perp = get_perpendicular_fluctuation(B0, dB)

print(f"Background B Field (B0): {np.array(B0)}")
print(f"Original Fluctuation (dB): {np.array(dB)}")
print("-" * 40)
print("To calculate helicity, one must use the components of the fluctuation")
print("that are perpendicular to the background B field.")
print(f"\nCalculated Perpendicular Fluctuation (dB_perp): {dB_perp}")
print("\nThe components of *this* vector (dB_perp) should be used to calculate helicity,")
print("NOT the original Y and Z components of the fluctuation.")
