import numpy as np

def simulate_resonance_measurement():
    """
    Simulates measuring a resonant coil with a VNA and loop probe,
    demonstrating the effect of coupling strength.
    """

    # --- Simulation Parameters ---
    # These represent the properties of the MRI coil and the measurement setup.
    f0 = 63.87  # Resonant frequency in MHz (e.g., for a 1.5T MRI)
    Q = 150.0   # Quality factor of the coil (how sharp the resonance is)
    measurement_noise_dB = 0.05 # A realistic noise floor for a VNA measurement in dB

    print("--- Simulation of MRI Coil Resonance Measurement ---\n")
    print(f"A functional coil is assumed to be resonant at f0 = {f0} MHz with a Quality Factor Q = {Q}.")
    print(f"Our measurement setup has a typical noise level of +/- {measurement_noise_dB} dB.\n")
    print("We will now simulate the measurement with two different coupling strengths.")

    # The S11 reflection dip (in dB) is modeled based on the coupling factor 'k'.
    # A larger 'k' means stronger coupling and a deeper, more visible dip.
    # The equation used is: S11_dB = 20 * log10(1 - k)
    # where k represents the fraction of power absorbed by the coil at resonance.

    # --- Case 1: Strong Coupling ---
    # This is what you see on coils where the resonance is visible.
    k_strong = 0.2
    s11_ideal_strong = 20 * np.log10(1 - k_strong)
    # Add random noise to simulate a real measurement
    s11_measured_strong = s11_ideal_strong + (np.random.rand() - 0.5) * 2 * measurement_noise_dB

    print("\n--- CASE 1: Strong Coupling ---")
    print("This occurs when the probe is close to an unshielded part of the coil.")
    print(f"Using a strong coupling factor k = {k_strong}:")
    print(f"The ideal S11 dip at {f0} MHz would be {s11_ideal_strong:.2f} dB.")
    print(f"A real measurement, including noise, might read around {s11_measured_strong:.2f} dB.")
    print("Result: This dip is very deep compared to the noise and is clearly visible.\n")


    # --- Case 2: Weak Coupling ---
    # This is the case for your shielded or geometrically complex coils.
    k_weak = 0.002
    s11_ideal_weak = 20 * np.log10(1 - k_weak)
    # Add random noise to simulate a real measurement
    s11_measured_weak = s11_ideal_weak + (np.random.rand() - 0.5) * 2 * measurement_noise_dB

    print("\n--- CASE 2: Weak Coupling ---")
    print("This occurs with shielded coils or when the probe is in a field null.")
    print(f"Using a weak coupling factor k = {k_weak}:")
    print(f"The ideal S11 dip at {f0} MHz would be only {s11_ideal_weak:.4f} dB.")
    print(f"A real measurement, including noise, might read around {s11_measured_weak:.4f} dB.")
    print(f"Result: This dip is much shallower than the measurement noise of +/- {measurement_noise_dB} dB and is therefore completely invisible.")
    print("\nThe resonance exists, but the coupling is too poor for the VNA to detect it.")


if __name__ == "__main__":
    simulate_resonance_measurement()