import numpy as np

def calculate_mfa_components():
    """
    Demonstrates the transformation of a magnetic field fluctuation vector from a
    GSE-like coordinate system to a Mean-Field-Aligned (MFA) system.
    
    This shows why using components perpendicular to the radial direction (By, Bz)
    is often an incorrect approximation for calculating helicity.
    """
    
    print("Step 1: Define the physical setup in GSE coordinates (units of nT).")
    # The X-axis points from the Sun to the Earth (radial).
    # The background Interplanetary Magnetic Field (IMF) at L1 is typically
    # not radial due to the Parker Spiral. We'll assume a ~45 degree angle
    # in the ecliptic (X-Y) plane.
    # B_avg = (Bx, By, Bz)
    B_avg = np.array([-4.0, 4.0, 0.5])
    
    # Let's define an instantaneous fluctuation vector, 'b'.
    # For an AIC wave, this fluctuation is mostly perpendicular to B_avg.
    b_fluctuation = np.array([0.5, 0.6, 2.0])

    print(f"Average Background Field B_avg (GSE): {B_avg} nT")
    print(f"Instantaneous Fluctuation b (GSE):   {b_fluctuation} nT\n")

    print("Step 2: Calculate the basis vectors for the Mean-Field-Aligned (MFA) frame.")
    # The parallel direction is the unit vector of the average magnetic field.
    e_parallel = B_avg / np.linalg.norm(B_avg)
    
    # To create the perpendicular plane, we use cross products.
    # We take the cross product of the parallel vector with the GSE Z-axis to get one perp vector.
    # This choice is arbitrary but standard; any non-collinear vector would work.
    e_perp2 = np.cross(e_parallel, np.array([0, 0, 1]))
    e_perp2 = e_perp2 / np.linalg.norm(e_perp2)
    
    # The third basis vector is found by the cross product of the other two to ensure a right-handed system.
    e_perp1 = np.cross(e_perp2, e_parallel)
    
    print(f"MFA Parallel unit vector (e_par):   {np.round(e_parallel, 3)}")
    print(f"MFA Perpendicular unit vector (e_p1): {np.round(e_perp1, 3)}")
    print(f"MFA Perpendicular unit vector (e_p2): {np.round(e_perp2, 3)}\n")

    print("Step 3: Project the fluctuation vector 'b' onto the new MFA basis.")
    b_p1 = np.dot(b_fluctuation, e_perp1)
    b_p2 = np.dot(b_fluctuation, e_perp2)
    b_par = np.dot(b_fluctuation, e_parallel)
    
    print("Fluctuation components in the correct MFA frame:")
    print(f"  b_perp1: {b_p1:.3f} nT")
    print(f"  b_perp2: {b_p2:.3f} nT")
    print(f"  b_parallel: {b_par:.3f} nT\n")
    
    print("Step 4: Compare perpendicular power calculation.")
    print("The denominator of the normalized helicity is the perpendicular power.")
    
    # Incorrect calculation using GSE components (assuming X is parallel)
    b_y_gse = b_fluctuation[1]
    b_z_gse = b_fluctuation[2]
    power_approx = b_y_gse**2 + b_z_gse**2
    
    print("INCORRECT Calculation using components perp to radial (GSE Y, Z):")
    print(f"  Perp Power = b_y^2 + b_z^2")
    print(f"  Perp Power = {b_y_gse:.3f}^2 + {b_z_gse:.3f}^2 = {power_approx:.3f} nT^2\n")

    # Correct calculation using the proper MFA perpendicular components
    power_correct = b_p1**2 + b_p2**2

    print("CORRECT Calculation using components perp to local B (MFA):")
    print(f"  Perp Power = b_perp1^2 + b_perp2^2")
    print(f"  Perp Power = {b_p1:.3f}^2 + {b_p2:.3f}^2 = {power_correct:.3f} nT^2")
    
    # Note: Total power is conserved in the transformation, as expected.
    # total_power_gse = np.sum(b_fluctuation**2)
    # total_power_mfa = b_p1**2 + b_p2**2 + b_par**2
    # print(f"\nDebug: Total power GSE={total_power_gse:.3f}, MFA={total_power_mfa:.3f}")

if __name__ == '__main__':
    calculate_mfa_components()
