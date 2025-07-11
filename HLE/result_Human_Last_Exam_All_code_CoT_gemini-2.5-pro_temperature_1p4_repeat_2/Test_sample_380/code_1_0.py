import numpy as np

def calculate_field_aligned_components():
    """
    Demonstrates the transformation to a Field-Aligned Coordinate (FAC) system
    for accurate helicity calculation.

    The user is correct that for waves propagating along the magnetic field,
    one must use components perpendicular to the local magnetic field, not
    just perpendicular to the radial direction. This script shows how to
    find those proper components.
    """
    # Step 1: Define sample magnetic field data in a standard RTN system (nT).
    # Bx is Radial, By is Tangential, Bz is Normal.
    # This sample represents a typical Parker Spiral field at L1, which is NOT radial.
    # We create a time series with some fluctuations around a mean value.
    np.random.seed(0)
    time_points = 1000
    # Mean field (B0) components
    b0_x, b0_y, b0_z = 3.0, -3.0, 0.5
    # Add some fluctuations (simulating waves)
    bx_series = b0_x + 0.1 * np.random.randn(time_points)
    by_series = b0_y + 0.5 * np.random.randn(time_points)
    bz_series = b0_z + 0.5 * np.random.randn(time_points)

    # Step 2: Calculate the mean magnetic field vector (B0) over the interval.
    # This vector defines the local background field direction.
    B0 = np.array([np.mean(bx_series), np.mean(by_series), np.mean(bz_series)])
    B0_magnitude = np.linalg.norm(B0)
    
    # Let's check if the field is radial. We can calculate the Parker Spiral Angle.
    # The angle is between the projection of B in the XY (ecliptic) plane and the radial (X) direction.
    parker_angle_rad = np.arctan2(B0[1], B0[0])
    parker_angle_deg = np.degrees(parker_angle_rad)

    print("--- Analysis of the Local Magnetic Field ---")
    print(f"Average Magnetic Field Vector B0 (RTN): {B0[0]:.2f}i {B0[1]:.2f}j + {B0[2]:.2f}k nT")
    print(f"Parker Spiral Angle: {parker_angle_deg:.2f} degrees")
    print("Conclusion: Since the angle is not 0, the field is NOT radial.\n")

    # Step 3: Construct the Field-Aligned Coordinate (FAC) system.
    # The system will have one axis parallel to B0, and two perpendicular axes.
    # e_parallel is the unit vector along B0.
    e_parallel = B0 / B0_magnitude

    # The first perpendicular axis can be defined using a cross product with a
    # stable direction, like the radial (X) axis [1, 0, 0].
    # This defines a plane perpendicular to B0.
    x_radial_unit_vector = np.array([1., 0., 0.])
    e_perp1 = np.cross(e_parallel, x_radial_unit_vector)
    e_perp1 = e_perp1 / np.linalg.norm(e_perp1)

    # The second perpendicular axis completes the right-handed system.
    e_perp2 = np.cross(e_parallel, e_perp1)
    
    print("--- Field-Aligned Coordinate (FAC) System ---")
    print("To correctly calculate helicity, we must use components perpendicular to B0.")
    print("These are the unit vectors for the two perpendicular directions:")
    print(f"Perpendicular Direction 1 (e_perp1): [{e_perp1[0]:.3f}, {e_perp1[1]:.3f}, {e_perp1[2]:.3f}]")
    print(f"Perpendicular Direction 2 (e_perp2): [{e_perp2[0]:.3f}, {e_perp2[1]:.3f}, {e_perp2[2]:.3f}]")
    print(f"Parallel Direction (e_parallel):      [{e_parallel[0]:.3f}, {e_parallel[1]:.3f}, {e_parallel[2]:.3f}]\n")

    # Step 4: Express the correct formula for normalized magnetic helicity.
    # The magnetic field fluctuations (delta_B) should be projected onto these
    # new perpendicular axes to get delta_b_perp1 and delta_b_perp2.
    # delta_b_perp1 = np.dot(delta_B, e_perp1)
    # delta_b_perp2 = np.dot(delta_B, e_perp2)
    # The helicity is then calculated from the power spectra of these new components.

    perp1_name = "B_perp1"
    perp2_name = "B_perp2"
    
    print("--- Correct Helicity Equation ---")
    print("The magnetic field fluctuations should be projected onto the new perpendicular axes.")
    print(f"The resulting fluctuating components are '{perp1_name}' and '{perp2_name}'.")
    print("The normalized magnetic helicity (sigma_m) is then defined as:")
    
    # Print the final equation using the new component names
    print(f"\n sigma_m(f) = (2 * Im[P_({perp1_name}{perp2_name})(f)]) / (P_({perp1_name}{perp1_name})(f) + P_({perp2_name}{perp2_name})(f))\n")
    
    print("Where:")
    print(f" P_({perp1_name}{perp2_name}) is the cross-power spectral density of the two perpendicular components.")
    print(f" P_({perp1_name}{perp1_name}) and P_({perp2_name}{perp2_name}) are the power spectral densities.")
    print(" Im[...] denotes the imaginary part.")

# Execute the function
calculate_field_aligned_components()

<<<Your physical intuition is correct. For a rigorous analysis of Alfven Ion Cyclotron waves, one must calculate the normalized magnetic helicity using magnetic field components that are perpendicular to the local background magnetic field, not just perpendicular to the radial direction. The local magnetic field at the L1 point is generally not radial due to the Parker Spiral. Therefore, a transformation into a Field-Aligned Coordinate (FAC) system is the proper method.>>>