def calculate_larmor_frequency(B0_values):
    """
    Calculates the Larmor frequency for hydrogen protons given scanner magnetic field strengths.

    The Larmor equation is f = γ * B₀, where:
    - f is the Larmor frequency in Hertz (Hz).
    - γ (gamma) is the gyromagnetic ratio, a constant for a given nucleus.
      For hydrogen protons (¹H), γ ≈ 42.5774 MHz/T.
    - B₀ is the external magnetic field strength in Tesla (T).
    """
    # Gyromagnetic ratio for Protons (¹H) in MHz per Tesla
    gamma_mhz_per_tesla = 42.5774

    print("Calculating Larmor frequency (f) for hydrogen protons using: f = γ * B₀")
    print(f"Using gyromagnetic ratio (γ) for ¹H = {gamma_mhz_per_tesla:.4f} MHz/T\n")

    for b0 in B0_values:
        # Calculate frequency in MHz
        frequency_mhz = gamma_mhz_per_tesla * b0
        print(f"For a scanner with B₀ = {b0:.1f} T:")
        # The final output prints the equation with all numbers as requested
        print(f"  The equation is: {gamma_mhz_per_tesla:.4f} MHz/T * {b0:.1f} T = {frequency_mhz:.2f} MHz")
        print("-" * 30)

if __name__ == "__main__":
    # Common clinical MRI scanner field strengths in Tesla
    common_field_strengths = [1.5, 3.0, 7.0]
    calculate_larmor_frequency(common_field_strengths)
