import numpy as np

def explain_helicity_coordinates():
    """
    Illustrates the different coordinate systems for analyzing magnetic field fluctuations
    in the solar wind at L1, providing justification for the common use of a
    radial-based system.
    """
    print("This script illustrates the coordinate systems for magnetic helicity calculations.\n")

    # --- 1. Define Background Conditions at L1 ---
    # The solar wind flows radially from the Sun. We define this as the X-direction.
    radial_direction = np.array([1., 0., 0.])
    print(f"Solar wind bulk flow (radial) direction: {radial_direction}\n")

    # The background Interplanetary Magnetic Field (IMF) is not radial due to the Parker Spiral.
    # It typically makes a ~45-degree angle with the radial direction in the ecliptic plane.
    parker_angle_deg = 45.0
    parker_angle_rad = np.deg2rad(parker_angle_deg)
    # Let's define the mean magnetic field B vector (in nT)
    B_magnitude = 5.0
    B_mean = B_magnitude * np.array([np.cos(parker_angle_rad), np.sin(parker_angle_rad), 0.])
    print(f"Assumed local mean magnetic field B vector (nT): {np.round(B_mean, 3)}")
    print("--> Key point: The local magnetic field is NOT radial.\n")

    # Let's consider a single sample of a magnetic field fluctuation (the wave component).
    dB = np.array([0.1, 0.5, -0.2])
    print(f"A sample magnetic field fluctuation vector dB (nT): {dB}\n")
    print("="*60)

    # --- 2. Method 1: Radial Coordinate System (Standard for statistical studies) ---
    print("Method 1: Using a Radial-based Coordinate System (e.g., GSE)\n")
    print("Justification: Taylor's 'Frozen-in Flow' Hypothesis. Because the solar wind")
    print("flow is fast and radial, a time series at a spacecraft is treated as a spatial")
    print("scan along the radial direction.")
    # In this system, the perpendicular components for helicity are simply the Y and Z components.
    dB_perp_Y_radial = dB[1]
    dB_perp_Z_radial = dB[2]
    print("\nThe components perpendicular to the RADIAL direction are:")
    print(f"  Component 1 (dBy) = {dB_perp_Y_radial}")
    print(f"  Component 2 (dBz) = {dB_perp_Z_radial}\n")
    print("="*60)

    # --- 3. Method 2: Local Magnetic Field-Aligned Coordinate System ---
    print("Method 2: Using a Local B-field-aligned Coordinate System\n")
    print("Justification: This is physically more accurate for a single wave event, as")
    print("AIC waves propagate along the local magnetic field B.")
    # To find the component of dB perpendicular to B_mean, we project dB onto B_mean.
    # First, get a unit vector in the direction of the mean magnetic field.
    B_mean_unit = B_mean / np.linalg.norm(B_mean)
    # The parallel fluctuation is the dot product of dB with the unit vector.
    dB_parallel_scalar = np.dot(dB, B_mean_unit)
    # The vector for the parallel component is this scalar times the unit vector.
    dB_parallel_vector = dB_parallel_scalar * B_mean_unit
    # The perpendicular fluctuation vector is the total fluctuation minus the parallel part.
    dB_perp_vector = dB - dB_parallel_vector
    print(f"The component of dB parallel to B has magnitude: {np.round(dB_parallel_scalar, 3)}")
    print(f"The vector of dB perpendicular to B is: {np.round(dB_perp_vector, 3)}")
    print("\nTo calculate helicity, one would then define two orthogonal axes in the plane")
    print("perpendicular to B. This is more complex and why Method 1 is often preferred.")
    print("="*60)

if __name__ == '__main__':
    explain_helicity_coordinates()