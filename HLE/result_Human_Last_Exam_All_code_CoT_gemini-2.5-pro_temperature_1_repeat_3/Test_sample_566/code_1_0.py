import numpy as np

def simulate_broadband_pump_cars():
    """
    Simulates broadband CARS with a broadband pump to demonstrate the generation
    of a distinguishable anti-Stokes spectrum.
    
    In CARS microscopy:
    - A pump beam (wp) and Stokes beam (ws) excite a molecular vibration (Omega).
    - Resonance occurs when Omega = wp - ws.
    - A probe beam (wpr) scatters off this vibration to create an anti-Stokes signal (was).
    - The anti-Stokes signal appears at the frequency: was = wpr + Omega.

    This simulation uses a broadband pump, meaning 'wp' is a range of frequencies.
    """

    # --- Define Simulation Parameters (in arbitrary frequency units, e.g., cm^-1) ---

    # 1. Define the vibrational frequencies of a hypothetical molecule sample.
    # These are the distinct pieces of information we want to detect.
    molecular_vibrations_Omega = [2850, 2920, 3010]

    # 2. Define the laser beam properties.
    # Stokes beam is narrowband (a single frequency).
    ws_stokes = 10000
    
    # Probe beam is narrowband (a single frequency).
    wpr_probe = 12900 
    
    # Pump beam is broadband (a continuous range of frequencies).
    wp_pump_range = np.arange(12800, 13101, 1) # A range from 12800 to 13100

    print("--- Simulation Setup ---")
    print(f"Molecular Vibrations to Detect (Omega): {molecular_vibrations_Omega}")
    print(f"Narrowband Stokes Beam (ws): {ws_stokes}")
    print(f"Narrowband Probe Beam (wpr): {wpr_probe}")
    print(f"Broadband Pump Beam (wp) covers frequencies from {wp_pump_range[0]} to {wp_pump_range[-1]}.\n")
    
    print("--- Finding Resonances and Generated Anti-Stokes Signals ---")
    generated_anti_stokes_signals = []

    # Iterate through each molecular vibration to see if it can be excited.
    for omega in molecular_vibrations_Omega:
        # For a given vibration (Omega), the required pump frequency (wp) to achieve
        # resonance (Omega = wp - ws) is:
        required_wp = omega + ws_stokes
        
        # Check if this required pump frequency exists within our broadband pump beam.
        if required_wp in wp_pump_range:
            # If yes, a coherent vibration is excited. The probe beam scatters off it.
            # Calculate the frequency of the generated anti-Stokes signal.
            was_anti_stokes = wpr_probe + omega
            generated_anti_stokes_signals.append(was_anti_stokes)
            
            print(f"Vibration Omega = {omega}:")
            print(f"  - Is excited by pump frequency component wp = {required_wp}.")
            print(f"  - This generates a distinguishable anti-Stokes signal at was = wpr + Omega = {wpr_probe} + {omega} = {was_anti_stokes}.")

    print("\n--- Final Result ---")
    if generated_anti_stokes_signals:
        print("The broadband pump successfully generated an anti-Stokes signal.")
        print("The full anti-Stokes signal contains the following distinguishable frequencies:")
        print(sorted(generated_anti_stokes_signals))
        print("\nThis demonstrates that a broadband pump can generate an anti-Stokes beam containing distinguishable vibrational information.")
    else:
        print("No anti-Stokes signals were generated with the current laser tuning.")

simulate_broadband_pump_cars()