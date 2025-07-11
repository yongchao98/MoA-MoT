import numpy as np

def demonstrate_coordinate_system_effect():
    """
    Demonstrates the importance of using a Field-Aligned Coordinate (FAC) system
    for analyzing magnetic field fluctuations.
    """
    # --- 1. Simulation Setup ---
    # Define parameters for a realistic solar wind scenario at L1
    B0_magnitude = 5.0  # Mean magnetic field magnitude in nT
    # Parker Spiral angle at 1 AU is typically ~45 degrees
    parker_angle_deg = 45.0
    parker_angle_rad = np.deg2rad(parker_angle_deg)

    # Wave properties (simulating an AIC wave)
    wave_amplitude = 1.0  # Wave amplitude in nT
    wave_frequency = 0.5  # Wave frequency in Hz (arbitrary for demonstration)

    # Time array for the simulation
    t = np.linspace(0, 10, 1000)

    # --- 2. Generate Simulated Magnetic Field Data ---
    # The data is generated in a Radial-Tangential-Normal (RTN) system,
    # where R (or X) is the radial direction.
    
    # Calculate mean background field components based on the Parker Spiral angle
    B_mean_R = B0_magnitude * np.cos(parker_angle_rad)
    B_mean_T = B0_magnitude * np.sin(parker_angle_rad)
    B_mean_N = 0.0

    # Generate transverse wave fluctuations. For this example, let's create a
    # circularly polarized wave in the T-N plane (perpendicular to radial).
    delta_B_T = wave_amplitude * np.cos(2 * np.pi * wave_frequency * t)
    delta_B_N = wave_amplitude * np.sin(2 * np.pi * wave_frequency * t)
    delta_B_R = np.zeros_like(t) # No radial fluctuation component

    # Total measured magnetic field (background + wave)
    B_R = B_mean_R + delta_B_R
    B_T = B_mean_T + delta_B_T
    B_N = B_mean_N + delta_B_N

    # --- 3. Analysis (The "Incorrect" Way) ---
    # Calculate fluctuating power perpendicular to the RADIAL direction.
    # This assumes the (T, N) components are the correct ones to use.
    power_perp_radial = np.mean(delta_B_T**2 + delta_B_N**2)

    # --- 4. Analysis (The "Correct" Way) ---
    # First, find the true mean magnetic field vector
    B_mean_vector = np.array([np.mean(B_R), np.mean(B_T), np.mean(B_N)])

    # Second, define the Field-Aligned Coordinate (FAC) system basis vectors
    # b_parallel is a unit vector along the mean magnetic field
    b_parallel = B_mean_vector / np.linalg.norm(B_mean_vector)

    # Create the perpendicular basis vectors. We ensure they are orthogonal.
    # A common method is to use a fixed vector (e.g., [0,0,1]) and cross products.
    b_perp2 = np.cross(b_parallel, [1, 0, 0]) # Cross with R-hat
    b_perp2 /= np.linalg.norm(b_perp2)
    b_perp1 = np.cross(b_perp2, b_parallel)
    
    # Third, project the FLUCTUATING part of the field onto the FAC system.
    # Fluctuating field vectors as a list of [dR, dT, dN]
    delta_B_vectors = np.vstack((delta_B_R, delta_B_T, delta_B_N)).T

    # Project onto perpendicular directions using dot product
    delta_B_perp1 = np.dot(delta_B_vectors, b_perp1)
    delta_B_perp2 = np.dot(delta_B_vectors, b_perp2)

    # Calculate fluctuating power perpendicular to the MEAN MAGNETIC FIELD
    power_perp_B0 = np.mean(delta_B_perp1**2 + delta_B_perp2**2)
    
    # We can also calculate the power of the compressional component (parallel fluctuation)
    delta_B_parallel = np.dot(delta_B_vectors, b_parallel)
    power_parallel_B0 = np.mean(delta_B_parallel**2)

    # --- 5. Print Results and Compare ---
    print("--- Simulation Parameters ---")
    print(f"Mean Magnetic Field Magnitude |B0| = {B0_magnitude:.2f} nT")
    print(f"Parker Spiral Angle = {parker_angle_deg:.1f} degrees")
    print(f"True Wave Power (Amplitude^2) = {wave_amplitude**2:.2f} nT^2")
    print("\n--- Mean Field Vector (RTN coordinates) ---")
    print(f"B_mean_R = {B_mean_vector[0]:.2f} nT")
    print(f"B_mean_T = {B_mean_vector[1]:.2f} nT")
    print(f"B_mean_N = {B_mean_vector[2]:.2f} nT")
    
    print("\n--- Power Calculation Comparison ---")
    print("The final two numbers show the results of the power calculation.")
    print("Magnetic helicity would be calculated from the fluctuations that constitute these powers.")
    print("\nMETHOD 1 (Incorrect Approximation): Power perpendicular to the RADIAL axis")
    # This power is calculated from delta_B_T**2 + delta_B_N**2
    print(f"Power_perp_radial = {power_perp_radial:.4f} nT^2")

    print("\nMETHOD 2 (Physically Correct): Power perpendicular to the MEAN MAGNETIC FIELD axis")
    # This power is calculated from delta_B_perp1**2 + delta_B_perp2**2
    print(f"Power_perp_B0 = {power_perp_B0:.4f} nT^2")
    
    print("\n--- Conclusion ---")
    print(f"Notice the calculated power differs. The 'correct' transverse power is {power_perp_B0:.2f} nT^2, not {power_perp_radial:.2f} nT^2.")
    print(f"The missing power ({power_parallel_B0:.2f} nT^2) now incorrectly appears as a compressional (parallel) fluctuation.")
    print("This demonstrates that for accurate physical quantities like helicity, a Field-Aligned Coordinate system must be used.")


if __name__ == '__main__':
    demonstrate_coordinate_system_effect()
