import numpy as np

def transform_to_mfa(br_series, bt_series, bn_series):
    """
    Transforms magnetic field time series from RTN coordinates to a
    Mean Field Aligned (MFA) coordinate system.

    Args:
        br_series (np.array): Time series of the radial magnetic field component.
        bt_series (np.array): Time series of the tangential magnetic field component.
        bn_series (np.array): Time series of the normal magnetic field component.

    Returns:
        tuple: A tuple containing the transformed magnetic field fluctuations:
               (b_perp1, b_perp2, b_parallel)
    """
    # 1. Calculate the mean magnetic field vector <B>
    # This vector defines the direction of the local magnetic field.
    mean_br = np.mean(br_series)
    mean_bt = np.mean(bt_series)
    mean_bn = np.mean(bn_series)
    B_mean_vec = np.array([mean_br, mean_bt, mean_bn])
    magnitude_B_mean = np.linalg.norm(B_mean_vec)

    print(f"Original Mean Field Vector <B> in RTN (nT):")
    print(f"  <Br> = {B_mean_vec[0]:.4f}")
    print(f"  <Bt> = {B_mean_vec[1]:.4f}")
    print(f"  <Bn> = {B_mean_vec[2]:.4f}\n")


    # 2. Define the basis vectors for the MFA coordinate system
    #    The parallel direction is aligned with the mean magnetic field.
    e_parallel = B_mean_vec / magnitude_B_mean

    #    To create the perpendicular axes, we need a vector that is not
    #    collinear with e_parallel. The radial vector is a good choice.
    R_vec = np.array([1.0, 0.0, 0.0])

    #    The first perpendicular direction is found using the cross product.
    e_perp2 = np.cross(e_parallel, R_vec)
    e_perp2 /= np.linalg.norm(e_perp2)

    #    The second perpendicular direction completes the right-handed system.
    e_perp1 = np.cross(e_perp2, e_parallel)

    print(f"MFA Coordinate System Basis Vectors (in RTN frame):")
    print(f"  e_parallel = [{e_parallel[0]:.4f}, {e_parallel[1]:.4f}, {e_parallel[2]:.4f}]")
    print(f"  e_perp1    = [{e_perp1[0]:.4f}, {e_perp1[1]:.4f}, {e_perp1[2]:.4f}]")
    print(f"  e_perp2    = [{e_perp2[0]:.4f}, {e_perp2[1]:.4f}, {e_perp2[2]:.4f}]\n")

    # 3. Calculate magnetic field fluctuations by subtracting the mean field
    b_fluctuations_r = br_series - mean_br
    b_fluctuations_t = bt_series - mean_bt
    b_fluctuations_n = bn_series - mean_bn

    # 4. Project the fluctuations onto the new MFA basis vectors
    #    These are the components to be used for helicity calculation.
    #    (b_perp1 and b_perp2)
    b_parallel = e_parallel[0] * b_fluctuations_r + e_parallel[1] * b_fluctuations_t + e_parallel[2] * b_fluctuations_n
    b_perp1 = e_perp1[0] * b_fluctuations_r + e_perp1[1] * b_fluctuations_t + e_perp1[2] * b_fluctuations_n
    b_perp2 = e_perp2[0] * b_fluctuations_r + e_perp2[1] * b_fluctuations_t + e_perp2[2] * b_fluctuations_n
    
    print("These are the physically correct components for calculating magnetic helicity:")
    print(" - b_perp1: Component in the first perpendicular direction.")
    print(" - b_perp2: Component in the second perpendicular direction.")

    return b_perp1, b_perp2, b_parallel

if __name__ == '__main__':
    # --- Example Usage ---
    # Create synthetic data simulating a non-radial magnetic field (Parker Spiral)
    # with a superimposed circularly polarized wave.
    time = np.linspace(0, 60, 600)  # 60 seconds of data

    # Mean field in RTN coordinates (e.g., [4 nT, 4 nT, 0.5 nT]), which is not radial.
    B0_r, B0_t, B0_n = 4.0, 4.0, 0.5

    # Wave fluctuations (circularly polarized in the T-N plane for simplicity)
    wave_amplitude = 0.5  # nT
    wave_frequency = 0.1  # Hz
    wave_bt = wave_amplitude * np.cos(2 * np.pi * wave_frequency * time)
    wave_bn = wave_amplitude * np.sin(2 * np.pi * wave_frequency * time) # left-hand polarized

    # Total measured field in RTN coordinates
    Br_total = B0_r + np.zeros_like(time)
    Bt_total = B0_t + wave_bt
    Bn_total = B0_n + wave_bn

    print("--- Transforming Magnetic Field Data to MFA Coordinates ---\n")
    # Perform the transformation
    b_perp1_ts, b_perp2_ts, b_parallel_ts = transform_to_mfa(Br_total, Bt_total, Bn_total)

    # You would now use b_perp1_ts and b_perp2_ts in the formula for
    # normalized magnetic helicity.
    print("\nExample transformation of the first data point:")
    B_first_point = np.array([Br_total[0], Bt_total[0], Bn_total[0]])
    b_fluctuation_first_point = B_first_point - np.array([np.mean(Br_total), np.mean(Bt_total), np.mean(Bn_total)])
    
    # We can reconstruct the MFA components for the first point
    B_mean_vec = np.array([np.mean(Br_total), np.mean(Bt_total), np.mean(Bn_total)])
    e_parallel_vec = B_mean_vec / np.linalg.norm(B_mean_vec)
    e_perp2_vec = np.cross(e_parallel_vec, np.array([1.,0.,0.]))
    e_perp2_vec /= np.linalg.norm(e_perp2_vec)
    e_perp1_vec = np.cross(e_perp2_vec, e_parallel_vec)
    
    bp1_first = np.dot(e_perp1_vec, b_fluctuation_first_point)
    bp2_first = np.dot(e_perp2_vec, b_fluctuation_first_point)
    bpar_first = np.dot(e_parallel_vec, b_fluctuation_first_point)
    
    print(f"Fluctuation vector at t=0 (RTN): [{b_fluctuation_first_point[0]:.4f}, {b_fluctuation_first_point[1]:.4f}, {b_fluctuation_first_point[2]:.4f}]")
    print(f"Projected onto MFA axes:         b_perp1={bp1_first:.4f}, b_perp2={bp2_first:.4f}, b_parallel={bpar_first:.4f}")
