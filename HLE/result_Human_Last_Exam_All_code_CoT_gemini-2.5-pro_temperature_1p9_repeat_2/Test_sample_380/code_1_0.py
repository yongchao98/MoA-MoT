import numpy as np
from scipy.signal import hilbert

def calculate_magnetic_helicity():
    """
    This function demonstrates the correct way to calculate normalized magnetic
    helicity by transforming magnetic field data into a Field-Aligned
    Coordinate (FAC) system.
    """
    # 1. Create synthetic magnetic field data in a radial-based system (RTN)
    # R-axis is radial from Sun, T is tangential, N is normal.
    # We will simulate a background field consistent with the Parker spiral at 1 AU,
    # where the field is at ~45 degrees to the radial direction.
    # Units are in nanoTesla (nT).
    B_mean_R = 3.5  # Radial component
    B_mean_T = -3.5 # Tangential component (creates a ~45 deg angle)
    B_mean_N = 0.5  # Small normal component
    
    # Create a synthetic left-hand polarized Alfven Ion Cyclotron (AIC) wave.
    # This wave propagates on top of the background field.
    sampling_rate = 10  # Hz
    duration = 100  # seconds
    num_points = int(duration * sampling_rate)
    time = np.linspace(0, duration, num_points)
    wave_frequency = 0.5  # Hz
    wave_amplitude = 0.5  # nT
    
    # Fluctuating components in the wave frame (b_perp1, b_perp2)
    # We define a left-hand polarized wave (b_perp2 leads b_perp1 by 90 deg)
    b_perp1 = wave_amplitude * np.cos(2 * np.pi * wave_frequency * time)
    b_perp2 = wave_amplitude * np.sin(2 * np.pi * wave_frequency * time)

    # Let's assume the wave propagates perfectly along the mean field direction.
    # The fluctuating field needs to be represented in the RTN system.
    # For this demonstration, we'll simplify and align the wave's perpendicular
    # plane roughly with the T-N plane for illustrative purposes. This is not
    # strictly necessary for the FAC transformation to work.
    b_R = np.zeros(num_points) # Wave is transverse, no fluctuation along propagation dir.
    b_T = b_perp1
    b_N = b_perp2

    # Total measured field in RTN coordinates
    B_R = B_mean_R + b_R
    B_T = B_mean_T + b_T
    B_N = B_mean_N + b_N

    print("--- Step 1: Input Data & Mean Field ---")
    B_mean_vec = np.array([np.mean(B_R), np.mean(B_T), np.mean(B_N)])
    print(f"Mean Magnetic Field Vector (RTN coords) B_mean = {B_mean_vec.round(2)} nT")
    print("Note: The field is not purely radial (first component is not the only non-zero one).\n")

    # 2. Define the Field-Aligned Coordinate (FAC) system
    # e_para: unit vector parallel to the mean magnetic field.
    e_para = B_mean_vec / np.linalg.norm(B_mean_vec)
    
    # e_perp2: unit vector perpendicular to the plane formed by e_para and the Radial direction.
    # This choice ensures one of the perpendicular axes lies in the R-T plane.
    e_perp2 = np.cross(e_para, np.array([1, 0, 0]))
    e_perp2 = e_perp2 / np.linalg.norm(e_perp2)

    # e_perp1: Completes the right-handed system.
    e_perp1 = np.cross(e_perp2, e_para)

    print("--- Step 2: Field-Aligned Coordinate (FAC) System Basis Vectors ---")
    print(f"Parallel axis (e_para):      {e_para.round(3)}")
    print(f"Perpendicular axis 1 (e_perp1): {e_perp1.round(3)}")
    print(f"Perpendicular axis 2 (e_perp2): {e_perp2.round(3)}\n")

    # 3. Transform the fluctuating magnetic field into the FAC system
    # First, get the fluctuating field by subtracting the mean
    b_fluctuations_RTN = np.vstack((B_R - B_mean_vec[0], B_T - B_mean_vec[1], B_N - B_mean_vec[2]))

    # Project the fluctuations onto the FAC basis vectors
    # b_fac_perp1/2 are the components needed for the helicity calculation
    b_fac_perp1 = np.dot(e_perp1, b_fluctuations_RTN)
    b_fac_perp2 = np.dot(e_perp2, b_fluctuations_RTN)
    
    print("--- Step 3: Calculate Normalized Magnetic Helicity (sigma_m) ---")
    # The normalized magnetic helicity, sigma_m, is defined as:
    # sigma_m = < (b_perp2 * H(b_perp1)) - (b_perp1 * H(b_perp2)) > / < b_perp1^2 + b_perp2^2 >
    # where H() is the Hilbert Transform and <...> is the time average.
    # A value of +1 is perfectly left-hand polarized; -1 is perfectly right-hand.

    # We use scipy.signal.hilbert. The analytic signal is `x + i*H(x)`.
    analytic_signal_1 = hilbert(b_fac_perp1)
    analytic_signal_2 = hilbert(b_fac_perp2)
    H_b1 = np.imag(analytic_signal_1)
    H_b2 = np.imag(analytic_signal_2)

    # Calculate the terms of the equation
    numerator_series = b_fac_perp2 * H_b1 - b_fac_perp1 * H_b2
    denominator_series = b_fac_perp1**2 + b_fac_perp2**2
    
    numerator = np.mean(numerator_series)
    denominator = np.mean(denominator_series)
    
    # Final calculation
    sigma_m = numerator / denominator

    print("The final equation is: sigma_m = Numerator / Denominator")
    print(f"Numerator = <b_perp2*H(b_perp1) - b_perp1*H(b_perp2)> = {numerator:.4f}")
    print(f"Denominator = <b_perp1^2 + b_perp2^2> = {denominator:.4f}")
    print(f"Normalized Magnetic Helicity (sigma_m) = {sigma_m:.4f}")
    print("\nResult is close to +1, confirming the left-hand polarization of the wave.")

if __name__ == '__main__':
    calculate_magnetic_helicity()
