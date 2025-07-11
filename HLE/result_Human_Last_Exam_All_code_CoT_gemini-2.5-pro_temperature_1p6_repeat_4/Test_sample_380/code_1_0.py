import numpy as np

def calculate_wave_power():
    """
    Demonstrates the difference between using a fixed radial coordinate system
    and a physically correct Field-Aligned Coordinate (FAC) system to analyze
    wave properties.
    """
    print("--- Setup: Creating a Synthetic Alfven Wave ---")

    # 1. Define the Mean Magnetic Field (B0) at a Parker Spiral angle of 45 degrees.
    # We use a Radial-Tangential-Normal (RTN) like system where X is Radial.
    # For a 45-degree angle, B_radial and B_tangential components have equal magnitude.
    # Let's say Bx = -5 nT and By = -5 nT (typical for 'away' magnetic sector).
    B0 = np.array([-5.0, -5.0, 0.0])
    B0_magnitude = np.linalg.norm(B0)
    print(f"Mean Magnetic Field B0 (Bx, By, Bz): {B0} nT")
    print(f"This field is at a 45-degree angle to the radial (X) direction.\n")

    # 2. Create a synthetic, left-hand circularly polarized wave.
    # The wave's true amplitude is 1 nT. Its power is amplitude^2 = 1.
    # This wave propagates along B0, so its fluctuations must be in the plane perpendicular to B0.
    
    # First, define the Field-Aligned Coordinate (FAC) basis vectors.
    e_parallel = B0 / B0_magnitude
    # Create a first perpendicular vector using a cross product with the Z-axis.
    e_perp1 = np.cross(e_parallel, np.array([0, 0, 1]))
    e_perp1 /= np.linalg.norm(e_perp1)
    # Create the second perpendicular vector to complete the right-handed system.
    e_perp2 = np.cross(e_parallel, e_perp1)

    print("--- Analysis: Calculating Perpendicular Wave Power ---")
    
    # Create wave signal over time 't'.
    t = np.linspace(0, 2 * np.pi, 100)
    wave_amplitude = 1.0
    
    # Fluctuating field 'b' in the perpendicular FAC plane
    b_perp1 = wave_amplitude * np.cos(t)
    b_perp2 = wave_amplitude * np.sin(t) # Circular polarization

    # The full fluctuating field vector in the original (X, Y, Z) system is:
    # b = b_perp1 * e_perp1 + b_perp2 * e_perp2
    b = np.outer(b_perp1, e_perp1) + np.outer(b_perp2, e_perp2)
    
    # The total measured magnetic field is the mean field plus the fluctuation.
    B_total = B0 + b

    # --- Method 1: The NAIVE approach (using components perp. to radial) ---
    # Here, we incorrectly use the fluctuations in the Y and Z directions.
    b_y_naive = B_total[:, 1] - B0[1]
    b_z_naive = B_total[:, 2] - B0[2]

    mean_sq_by = np.mean(b_y_naive**2)
    mean_sq_bz = np.mean(b_z_naive**2)
    power_naive = mean_sq_by + mean_sq_bz

    print("\nMETHOD 1: NAIVE (using Y and Z components of the original system)")
    print("This assumes the wave propagates along the X (radial) axis, which is incorrect.")
    print(f"Calculated Perpendicular Power = <b_y^2> + <b_z^2>")
    print(f"Result: {mean_sq_by:.2f} + {mean_sq_bz:.2f} = {power_naive:.2f} nT^2")
    print("This result incorrectly underestimates the true power.")


    # --- Method 2: The CORRECT approach (using FAC system) ---
    # Here, we use the perpendicular components we defined earlier.
    # The calculation is simply the mean square of the known perpendicular components.
    mean_sq_perp1 = np.mean(b_perp1**2)
    mean_sq_perp2 = np.mean(b_perp2**2)
    power_correct = mean_sq_perp1 + mean_sq_perp2
    
    print("\nMETHOD 2: CORRECT (using components in the Field-Aligned Coordinate system)")
    print("This uses components truly perpendicular to the mean magnetic field B0.")
    print(f"True Perpendicular Power = <b_perp1^2> + <b_perp2^2>")
    print(f"Result: {mean_sq_perp1:.2f} + {mean_sq_perp2:.2f} = {power_correct:.2f} nT^2")
    print("This result recovers the true wave power (Amplitude^2 = 1.0^2 = 1.0).")


if __name__ == '__main__':
    calculate_wave_power()