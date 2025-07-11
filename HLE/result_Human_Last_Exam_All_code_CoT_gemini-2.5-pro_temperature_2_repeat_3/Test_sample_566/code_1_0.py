import sys

def explain_broadband_cars():
    """
    This function explains the principle of broadband CARS microscopy and demonstrates
    why it produces a signal with distinguishable vibrational information.
    """
    print("Step 1: Define the laser beam frequencies for a CARS experiment.")
    print("Broadband CARS aims to excite multiple vibrations at once. We'll model this with a narrowband pump laser and a broadband Stokes laser which contains many frequencies.")
    print("Let's use typical frequencies in wavenumbers (cm^-1).\n")

    # Define a narrowband pump frequency
    omega_p = 12500  # Corresponds to an 800 nm laser

    # Define two target molecular vibrations (e.g., C-H stretches)
    Omega_1 = 2900
    Omega_2 = 3050

    # The required Stokes frequencies to excite these vibrations are:
    omega_S1 = omega_p - Omega_1
    omega_S2 = omega_p - Omega_2

    print(f"Pump Frequency (ω_p): {omega_p} cm^-1")
    print(f"Target Vibration 1 (Ω_1): {Omega_1} cm^-1")
    print(f"Target Vibration 2 (Ω_2): {Omega_2} cm^-1")
    print(f"Corresponding Stokes Frequency 1 (ω_S1): {omega_S1} cm^-1")
    print(f"Corresponding Stokes Frequency 2 (ω_S2): {omega_S2} cm^-1\n")

    print("Step 2: Calculate the Anti-Stokes signal frequencies.")
    print("The anti-Stokes signal (ω_as) is generated at ω_as = ω_p + Ω, where Ω = ω_p - ω_S.")
    print("Because we excite two different vibrations (Ω_1 and Ω_2), we will generate two distinct anti-Stokes frequencies.\n")

    # Calculate the resulting anti-Stokes frequencies
    omega_as1 = omega_p + Omega_1
    omega_as2 = omega_p + Omega_2

    print("--- Calculating the anti-Stokes signal for Vibration 1 ---")
    print(f"Final Equation: {omega_as1} = {omega_p} + {Omega_1}")
    print("Outputting each number in the final equation:")
    print(omega_p)
    print(Omega_1)
    print(omega_as1)
    print("----------------------------------------------------------\n")

    print("--- Calculating the anti-Stokes signal for Vibration 2 ---")
    print(f"Final Equation: {omega_as2} = {omega_p} + {Omega_2}")
    print("Outputting each number in the final equation:")
    print(omega_p)
    print(Omega_2)
    print(omega_as2)
    print("----------------------------------------------------------\n")
    
    print("Step 3: Conclusion")
    print(f"The resulting anti-Stokes signal contains distinct frequencies at {omega_as1} cm^-1 and {omega_as2} cm^-1.")
    print("These can be easily separated by a spectrometer, meaning the signal contains 'distinguishable information' about different molecules.")
    print("\nTherefore, doing broadband CARS generates an anti-Stokes beam that contains distinguishable vibrational information.")

if __name__ == '__main__':
    explain_broadband_cars()

<<<C>>>