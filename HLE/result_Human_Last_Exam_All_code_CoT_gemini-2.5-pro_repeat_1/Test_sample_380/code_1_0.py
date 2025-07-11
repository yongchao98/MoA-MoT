import numpy as np

def solve_helicity_components():
    """
    This script demonstrates the difference between using a simplified radial
    coordinate system and a more rigorous field-aligned coordinate system
    for calculating wave properties like magnetic helicity.
    """

    # --- Step 1: Define a realistic magnetic field at L1 ---
    # The Interplanetary Magnetic Field (IMF) at L1 is not purely radial.
    # Due to the Sun's rotation, it forms a Parker Spiral. At 1 AU (Earth's orbit),
    # the angle is typically ~45 degrees from the radial direction in the ecliptic plane.
    # We use a Geocentric Solar Ecliptic (GSE) like frame where X is Sun->Earth.
    # A negative Bx means the field points towards the Sun.
    # Let's define a mean magnetic field vector B_mean in nanoTeslas (nT).
    Bx_mean = -4.0  # Radial component
    By_mean = 4.0   # Azimuthal component (in the ecliptic plane)
    Bz_mean = 1.0   # North-South component (out of the ecliptic)
    B_mean = np.array([Bx_mean, By_mean, Bz_mean])

    # Now, let's add a small fluctuation (representing the AIC wave) to this mean field.
    # This is the instantaneous total magnetic field B_total that would be measured.
    dBx = 0.2
    dBy = 0.8
    dBz = -0.6
    B_fluctuation = np.array([dBx, dBy, dBz])
    B_total = B_mean + B_fluctuation

    print("--- Analysis of Magnetic Field Components for Helicity Calculation ---")
    print(f"The total measured magnetic field vector is B = <{B_total[0]:.2f}, {B_total[1]:.2f}, {B_total[2]:.2f}> nT.")
    print(f"The local mean magnetic field vector is B_0 = <{B_mean[0]:.2f}, {B_mean[1]:.2f}, {B_mean[2]:.2f}> nT.\n")


    # --- Step 2: Case 1 - The Simplified "Radial" Assumption ---
    # This method assumes the background field is purely radial (along the X-axis).
    # The helicity is then calculated using the components perpendicular to the X-axis,
    # which are simply the Y and Z components of the total measured field.
    B_perp_radial_approx = np.array([0, B_total[1], B_total[2]])

    print("--- Case 1: Simplified Radial Approximation ---")
    print("This assumes the main field is along the radial (X) axis.")
    print("The perpendicular field components used for helicity are taken from the Y-Z plane.")
    print("Equation for the perpendicular field vector B_perp_radial:")
    # The final equation as requested
    print(f"B_perp_radial = <0.00, {B_perp_radial_approx[1]:.2f}, {B_perp_radial_approx[2]:.2f}> nT\n")


    # --- Step 3: Case 2 - The Rigorous "Field-Aligned" Method ---
    # This is the physically correct method for field-aligned waves like AIC waves.
    # We must find the components of the wave's magnetic field (B_fluctuation) that are
    # truly perpendicular to the local background magnetic field (B_mean).

    # First, find the unit vector for the mean magnetic field direction.
    b_hat = B_mean / np.linalg.norm(B_mean)

    # Project the fluctuation vector onto the mean field direction to find its parallel component.
    B_fluc_parallel_scalar = np.dot(B_fluctuation, b_hat)
    B_fluc_parallel_vector = B_fluc_parallel_scalar * b_hat

    # The perpendicular component of the fluctuation is the total fluctuation minus its parallel part.
    B_fluc_perp_local = B_fluctuation - B_fluc_parallel_vector

    print("--- Case 2: Rigorous Field-Aligned Calculation ---")
    print("This correctly uses the components perpendicular to the local mean magnetic field B_0.")
    print("The perpendicular field components are found by projecting the wave fluctuation and taking the remainder.")
    print("Equation for the perpendicular field vector B_perp_local:")
    # The final equation as requested
    print(f"B_perp_local = <{B_fluc_perp_local[0]:.2f}, {B_fluc_perp_local[1]:.2f}, {B_fluc_perp_local[2]:.2f}> nT\n")

    print("--- Conclusion ---")
    print("As you can see, the vectors used for the helicity calculation are numerically different.")
    print("Using the Y and Z components is a simplification based on the measurement coordinate system (like GSE).")
    print("For rigorous analysis of field-aligned waves, one must transform to a field-aligned coordinate system.")


if __name__ == '__main__':
    solve_helicity_components()