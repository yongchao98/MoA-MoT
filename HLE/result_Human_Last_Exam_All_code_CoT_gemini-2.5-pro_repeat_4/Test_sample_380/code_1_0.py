import numpy as np

def explain_and_demonstrate_transformation():
    """
    Demonstrates the necessity of transforming to a Field-Aligned Coordinate (FAC)
    system to analyze waves propagating along the magnetic field.
    """
    print("--- Conceptual Demonstration of Coordinate Transformation for Wave Analysis ---")

    # Step 1: Define the environment in the measurement (GSE-like) frame.
    # The X-axis is assumed to be radial from the Sun.
    print("\nStep 1: Define the physical environment at the L1 point.")
    
    # B0 is the background magnetic field. We model it as a Parker Spiral at 45 degrees
    # in the XY-plane. It is NOT radial. A purely radial field would be [B, 0, 0].
    B0 = np.array([1.0, 1.0, 0.2])
    B0_magnitude = np.linalg.norm(B0)
    b_hat = B0 / B0_magnitude # Unit vector in the direction of B0

    # dB is the magnetic field of an AIC wave. By definition, it's a transverse wave,
    # so its magnetic field must be perpendicular to its direction of propagation (which is B0).
    # We can construct a vector perpendicular to B0 to represent the wave.
    # A simple way is to take the cross product of B0 with any non-parallel vector.
    # We'll use the Z-axis vector [0,0,1].
    wave_perp_vec1 = np.cross(b_hat, [0, 0, 1])
    wave_perp_vec1 /= np.linalg.norm(wave_perp_vec1) # Normalize
    
    # This represents an instantaneous snapshot of the wave's magnetic field.
    dB = wave_perp_vec1 
    
    print(f"Assumed Radial Direction (X-axis): [1, 0, 0]")
    print(f"Local Magnetic Field Direction (B_hat): {np.round(b_hat, 3)}")
    print(f"Instantaneous Wave Field (dB): {np.round(dB, 3)}")
    print(f"Dot product of dB and B_hat: {np.round(np.dot(dB, b_hat), 5)} (This is zero, confirming the wave is transverse to B0)")

    # Step 2: Analyze the wave in the original, incorrect frame.
    # If we incorrectly use components perpendicular to the *radial* direction (Y and Z),
    # we miss part of the wave signal.
    print("\nStep 2: Analyze the wave using components perpendicular to the RADIAL direction (Incorrect).")
    dB_x, dB_y, dB_z = dB
    print(f"The wave component along the radial (X) direction is dB_x = {dB_x:.3f}")
    print("Since dB_x is not zero, simply using dB_y and dB_z to calculate helicity would be incorrect and would miss this component.")

    # Step 3: Perform the correct coordinate transformation to a Field-Aligned Coordinate (FAC) system.
    print("\nStep 3: Transform to a Field-Aligned Coordinate (FAC) system (Correct).")
    
    # The parallel direction in FAC is the direction of the background magnetic field.
    z_fac = b_hat
    
    # We define the perpendicular plane. One axis can be the cross product of the B0 direction
    # and the radial direction (X-axis). This defines a perpendicular direction in the ecliptic plane.
    y_fac = np.cross(z_fac, [1, 0, 0])
    y_fac /= np.linalg.norm(y_fac)

    # The third axis completes the right-handed system.
    x_fac = np.cross(y_fac, z_fac)

    print(f"FAC Parallel Axis (z_fac, aligned with B0): {np.round(z_fac, 3)}")
    print(f"FAC Perpendicular Axis 1 (x_fac): {np.round(x_fac, 3)}")
    print(f"FAC Perpendicular Axis 2 (y_fac): {np.round(y_fac, 3)}")

    # Step 4: Express the wave vector in the new FAC basis and calculate helicity components.
    print("\nStep 4: Express the wave in FAC and identify components for helicity calculation.")
    
    # Project the wave vector dB onto the new FAC basis vectors.
    dB_fac_x = np.dot(dB, x_fac)
    dB_fac_y = np.dot(dB, y_fac)
    dB_fac_z = np.dot(dB, z_fac)
    
    print(f"Wave component in FAC parallel direction (dB_fac_z): {np.round(dB_fac_z, 5)}")
    print(f"Wave component in FAC perpendicular direction 1 (dB_fac_x): {np.round(dB_fac_x, 3)}")
    print(f"Wave component in FAC perpendicular direction 2 (dB_fac_y): {np.round(dB_fac_y, 3)}")
    print("\nConclusion: The wave field is entirely contained by the two perpendicular FAC components (dB_fac_x and dB_fac_y).")
    print("These are the components that must be used to correctly calculate the normalized magnetic helicity.")

if __name__ == '__main__':
    explain_and_demonstrate_transformation()