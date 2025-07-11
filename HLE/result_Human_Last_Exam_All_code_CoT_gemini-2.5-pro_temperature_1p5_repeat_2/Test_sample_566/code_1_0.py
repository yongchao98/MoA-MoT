import sys

def explain_broadband_cars():
    """
    Explains the principle of broadband CARS microscopy to answer the user's question.
    """
    print("### Analyzing Broadband CARS Microscopy ###\n")

    # Step 1: Explain the basic CARS process.
    print("1. The CARS (Coherent Anti-Stokes Raman Scattering) process uses multiple laser pulses to generate a signal.")
    print("   - It is enhanced when the frequency difference between a pump (ωp) and Stokes (ωs) beam matches a molecular vibration (ΩR).")
    print("   - A probe beam (usually the same as the pump) interacts with this vibration to create a signal.")

    # Step 2: Define the signal frequency and the core equation.
    print("\n2. The generated signal is at the anti-Stokes frequency (ω_as), given by the equation:")
    pump_freq_THz = 375  # Corresponds to 800 nm laser
    stokes_freq_THz = 288 # Corresponds to 1040 nm laser
    vibrational_freq_THz = pump_freq_THz - stokes_freq_THz
    anti_stokes_freq_THz = pump_freq_THz + vibrational_freq_THz
    print(f"   ω_as = ω_p + (ω_p - ω_s) = 2 * ω_p - ω_s")
    print(f"   For example: If pump is {pump_freq_THz} THz and Stokes is {stokes_freq_THz} THz...")
    print(f"   The vibrational frequency is {pump_freq_THz} - {stokes_freq_THz} = {vibrational_freq_THz} THz.")
    print(f"   The anti-Stokes signal is at {pump_freq_THz} + {vibrational_freq_THz} = {anti_stokes_freq_THz} THz.\n")

    # Step 3: Explain the "Broadband" aspect.
    print("3. 'Broadband' CARS aims to capture a full vibrational spectrum at once, unlike narrowband CARS which measures one vibration at a time.")
    print("   - This is achieved by using one laser with a very broad range of frequencies (a broadband beam) and one with a narrow range.")
    print("   - This excites a wide range of molecular vibrations simultaneously.\n")

    # Step 4: Describe the resulting signal and conclusion.
    print("4. This simultaneous excitation produces a broadband anti-Stokes signal (ω_as).")
    print("   - This signal is essentially a bundle of many different light frequencies, with each frequency corresponding to a specific molecular vibration.")
    print("   - By using a spectrometer to separate these frequencies, we can get a full vibrational spectrum.")
    print("   - Therefore, the anti-Stokes beam contains distinguishable information about the different chemical bonds in the sample.\n")

    # Step 5: Final conclusion based on the choices.
    print("### Conclusion ###")
    print("Based on this, the correct statement is that broadband CARS generates an anti-Stokes beam containing distinguishable information.")


if __name__ == '__main__':
    explain_broadband_cars()