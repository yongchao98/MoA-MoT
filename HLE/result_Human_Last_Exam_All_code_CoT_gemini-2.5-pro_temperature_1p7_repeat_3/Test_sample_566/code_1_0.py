import numpy as np

def simulate_broadband_pump_cars():
    """
    This script demonstrates how using a broadband pump beam in CARS microscopy
    can generate an anti-Stokes signal containing distinguishable vibrational information.
    """
    # --- Setup ---
    # We use wavenumbers (cm^-1) for all frequencies. Let's assume we have
    # a narrowband Stokes laser and a narrowband probe laser.
    omega_s = 10000.0  # Narrowband Stokes frequency in cm^-1 (~1000 nm)
    omega_probe = 12500.0 # Narrowband probe frequency in cm^-1 (~800 nm)

    # Let's assume our sample has several distinct molecular vibrational frequencies.
    # These are the pieces of information we want to distinguish in our final signal.
    # We will use some example Raman peaks for a common chemical like Toluene.
    Omega_vibrations = np.array([1004.0, 1211.0, 2917.0])

    print("--- Simulating Broadband Pump CARS ---")
    print(f"The system has a narrowband Stokes beam at ω_s = {omega_s} cm⁻¹.")
    print(f"It also has a narrowband probe beam at ω_probe = {omega_probe} cm⁻¹.")
    print(f"The pump beam is broadband, covering a wide range of frequencies.")
    print(f"The sample has the following vibrational modes (Ω_vib): {Omega_vibrations} cm⁻¹.\n")

    print("--- Calculation ---")
    print("The final anti-Stokes frequency (ω_as) is generated via the equation:")
    print("ω_as = ω_p_resonant + ω_probe - ω_s")
    print("where ω_p_resonant is the frequency from the broadband pump that excites a specific vibration.\n")

    # --- Main Calculation Loop ---
    # For each vibrational mode, we calculate the resulting anti-Stokes frequency.
    for i, Omega_vib in enumerate(Omega_vibrations):
        # From the broadband pump, the frequency that becomes resonant is ω_p_resonant = Ω_vib + ω_s
        omega_p_resonant = Omega_vib + omega_s

        # This resonant pump frequency generates the final signal
        omega_as = omega_p_resonant + omega_probe - omega_s

        print(f"For vibrational mode Ω_vib = {Omega_vib} cm⁻¹:")
        print(f"  1. The resonant pump frequency is: {Omega_vib} + {omega_s} = {omega_p_resonant} cm⁻¹")
        print(f"  2. The generated anti-Stokes signal frequency is:")
        print(f"     ω_as = {omega_p_resonant} + {omega_probe} - {omega_s} = {omega_as} cm⁻¹")
        print("-" * 40)

    print("\n--- Conclusion ---")
    print("Each unique vibrational frequency (Ω_vib) produces a unique and distinguishable")
    print("anti-Stokes frequency (ω_as). A spectrometer can measure this output spectrum,")
    print("revealing the sample's vibrational information all at once.")
    print("This confirms that you can generate an anti-Stokes beam that contains distinguishable information.")

# Run the simulation
simulate_broadband_pump_cars()