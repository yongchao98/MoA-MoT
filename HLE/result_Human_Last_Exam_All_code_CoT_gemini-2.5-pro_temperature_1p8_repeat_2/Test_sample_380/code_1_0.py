import numpy as np

def explain_and_demonstrate_coordinates():
    """
    Demonstrates the difference between analyzing wave components in a standard
    coordinate system vs. a field-aligned coordinate system.
    """
    # --- 1. Setup the Scenario ---
    # Define a background Interplanetary Magnetic Field (IMF) at the L1 point.
    # We use Geocentric Solar Ecliptic (GSE) coordinates, where X is from Earth to Sun.
    # The solar wind flows radially along the -X direction.
    # We model a Parker Spiral field with an angle of 45 degrees to the radial (X) axis.
    # Let's assume it lies in the X-Y plane for simplicity.
    parker_angle_deg = 45.0
    parker_angle_rad = np.deg2rad(parker_angle_deg)
    B_mag = 5.0  # Magnitude in nT
    
    # B_background is NOT radial. A purely radial field would be [Bx, 0, 0].
    B_background = np.array([B_mag * np.cos(parker_angle_rad), B_mag * np.sin(parker_angle_rad), 0.0])

    # Now, add a magnetic perturbation (like a part of a wave)
    # Let's assume a simple fluctuation purely in the Z-direction for clarity.
    b_fluctuation = np.array([0.0, 0.0, 1.0]) # 1 nT fluctuation

    # The total field measured by the spacecraft
    B_total = B_background + b_fluctuation
    
    print("--- Scenario Setup ---")
    print(f"Background IMF B_background = {np.round(B_background, 2)} nT")
    print(f"Fluctuation b_fluctuation = {np.round(b_fluctuation, 2)} nT")
    print(f"Total measured field B_total = {np.round(B_total, 2)} nT")
    print("-" * 30 + "\n")

    # --- 2. Common Method (Using GSE Y and Z components) ---
    # This method approximates the transverse components as being perpendicular
    # to the radial (X) direction.
    B_transverse_gse_y = B_total[1]
    B_transverse_gse_z = B_total[2]
    transverse_power_gse = B_transverse_gse_y**2 + B_transverse_gse_z**2
    
    print("--- Method 1: Common Approximation (GSE Frame) ---")
    print("Assumes transverse direction is perpendicular to the radial X-axis.")
    print(f"Component By = {B_transverse_gse_y:.2f} nT")
    print(f"Component Bz = {B_transverse_gse_z:.2f} nT")
    print(f"Apparent 'Transverse' Power (By^2 + Bz^2) = {transverse_power_gse:.2f} nT^2")
    print("-" * 30 + "\n")
    
    # --- 3. Rigorous Method (Using Field-Aligned Coordinates) ---
    # This method correctly identifies the components perpendicular to the LOCAL background field.
    # a. Define the Field-Aligned Coordinate (FAC) system basis vectors
    e_parallel = B_background / np.linalg.norm(B_background)
    
    # Create a second vector for the perpendicular plane, e.g., in the Z direction
    # if B_background is in XY plane.
    if np.abs(e_parallel[2]) < 0.9:
      e_perp1 = np.cross(e_parallel, [0, 0, 1])
      e_perp1 /= np.linalg.norm(e_perp1)
    else: # if B_background is aligned with Z, use Y
      e_perp1 = np.cross(e_parallel, [0, 1, 0])
      e_perp1 /= np.linalg.norm(e_perp1)
    
    e_perp2 = np.cross(e_parallel, e_perp1)

    # b. Project the FLUCTUATION onto the perpendicular directions
    b_perp1 = np.dot(b_fluctuation, e_perp1)
    b_perp2 = np.dot(b_fluctuation, e_perp2)
    b_parallel = np.dot(b_fluctuation, e_parallel)
    
    true_transverse_power = b_perp1**2 + b_perp2**2
    
    print("--- Method 2: Rigorous Physical Method (Field-Aligned Frame) ---")
    print("Identifies components truly perpendicular to the background B-field.")
    print(f"Unit vector parallel to B: e_|| = {np.round(e_parallel, 2)}")
    print(f"Perpendicular component 1 (b_perp1) = {b_perp1:.2f} nT")
    print(f"Perpendicular component 2 (b_perp2) = {b_perp2:.2f} nT")
    print(f"True Transverse Power (b_perp1^2 + b_perp2^2) = {true_transverse_power:.2f} nT^2")
    print("\n--- Conclusion ---")
    print("The power calculated by the two methods is different.")
    print("The GSE method (Method 1) incorrectly mixes parallel and perpendicular components.")
    print("For our example fluctuation purely in Z, the GSE 'transverse' power is 1.00 nT^2.")
    print(f"However, the true transverse power is only {true_transverse_power:.2f} nT^2 because part of the Z-fluctuation is parallel to the tilted background B-field.")
    
# Run the demonstration
explain_and_demonstrate_coordinates()