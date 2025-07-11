import numpy as np

def solve_helicity_coordinates():
    """
    Demonstrates the coordinate transformation required for helicity calculations.

    This script shows how to transform a magnetic field fluctuation vector from a
    standard RTN (Radial-Tangential-Normal) coordinate system to a Field-Aligned
    Coordinate (FAC) system, which is necessary for the physical analysis of
    field-aligned waves like AIC waves.
    """

    # --- Step 1: Define the physical setup in the RTN system ---
    # The RTN system is Sun-centered: R is Radial, T is Tangential, N is Normal.
    
    # Assume a typical average magnetic field (B0) at L1 (1 AU).
    # Total field strength is 5 nT, with a Parker Spiral angle of 45 degrees.
    # The field points towards the Sun along the spiral.
    b0_magnitude = 5.0  # nT
    parker_angle_deg = 45.0
    parker_angle_rad = np.deg2rad(parker_angle_deg)

    # B0 components in RTN (R points away from Sun, so B_R is negative)
    B0_R = -b0_magnitude * np.cos(parker_angle_rad)
    B0_T =  b0_magnitude * np.sin(parker_angle_rad)
    B0_N = 0.0 # Idealized in the ecliptic plane
    B0_rtn = np.array([B0_R, B0_T, B0_N])

    # Assume we measure a magnetic field fluctuation (delta_B) at a point in time.
    delta_B_rtn = np.array([0.5, -0.2, 0.8]) # An arbitrary fluctuation in nT

    print("--- Input Vectors in Radial-Tangential-Normal (RTN) System ---")
    print(f"Average Magnetic Field B0 (RTN): [{B0_rtn[0]:.4f}, {B0_rtn[1]:.4f}, {B0_rtn[2]:.4f}] nT")
    print(f"Fluctuation delta_B (RTN):       [{delta_B_rtn[0]:.4f}, {delta_B_rtn[1]:.4f}, {delta_B_rtn[2]:.4f}] nT\n")

    # --- Step 2: Create the basis vectors for the Field-Aligned Coordinate (FAC) system ---
    # The new system will have axes (perp1, perp2, parallel).
    
    # The 'parallel' axis is along the direction of the average magnetic field B0.
    e_parallel = B0_rtn / np.linalg.norm(B0_rtn)

    # To define the perpendicular plane, we use the RTN radial unit vector as a reference.
    r_hat = np.array([1.0, 0.0, 0.0])
    
    # The 'perp2' axis is perpendicular to both B0 and the radial direction.
    e_perp2 = np.cross(e_parallel, r_hat)
    e_perp2 = e_perp2 / np.linalg.norm(e_perp2)

    # The 'perp1' axis completes the right-handed system.
    e_perp1 = np.cross(e_perp2, e_parallel)

    print("--- Field-Aligned Coordinate (FAC) System Basis Vectors ---")
    print(f"e_parallel: [{e_parallel[0]:.4f}, {e_parallel[1]:.4f}, {e_parallel[2]:.4f}]")
    print(f"e_perp1:    [{e_perp1[0]:.4f}, {e_perp1[1]:.4f}, {e_perp1[2]:.4f}]")
    print(f"e_perp2:    [{e_perp2[0]:.4f}, {e_perp2[1]:.4f}, {e_perp2[2]:.4f}]\n")

    # --- Step 3: Project the fluctuation vector onto the new FAC basis ---
    # This is the coordinate transformation.
    delta_B_parallel = np.dot(delta_B_rtn, e_parallel)
    delta_B_perp1 = np.dot(delta_B_rtn, e_perp1)
    delta_B_perp2 = np.dot(delta_B_rtn, e_perp2)

    delta_B_fac = np.array([delta_B_perp1, delta_B_perp2, delta_B_parallel])

    print("--- Final Result: Fluctuation Vector in FAC System ---")
    print("The transformation is done by projecting delta_B onto the new basis vectors.")
    print("Equation: delta_B_fac = [dot(delta_B, e_perp1), dot(delta_B, e_perp2), dot(delta_B, e_parallel)]")
    print(f"Transformed delta_B (FAC): [{delta_B_fac[0]:.4f} (perp1), {delta_B_fac[1]:.4f} (perp2), {delta_B_fac[2]:.4f} (parallel)] nT\n")
    
    print("--- Conclusion ---")
    print("The normalized magnetic helicity for AIC waves would be calculated using the time series of")
    print("the perpendicular components: delta_B_perp1 and delta_B_perp2.")

solve_helicity_coordinates()