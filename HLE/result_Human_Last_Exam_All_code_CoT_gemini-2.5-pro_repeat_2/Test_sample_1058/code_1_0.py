# --- Experiment 1: Behavioral Data ---
rats_initial_study = 27
percentage_preferring_alcohol_initial = 15
percentage_preferring_sucrose_initial = 85

rats_larger_study = 967
percentage_preferring_alcohol_larger = 15.2

# --- Experiment 2: Electrophysiology Baseline ---
stimulation_intensity_uA = 40
ps_amplitude_alcohol_preferring_mV = -0.38
ps_amplitude_sucrose_preferring_mV = -0.17

# --- Experiment 3: Gene Knockdown (Slc6a11) ---
ps_amplitude_shRNA_knockdown_mV = -0.37
ps_amplitude_control_vector_mV = -0.16

# --- Output the key comparisons ---
print("--- Key Experimental Values ---")
print(f"Comparison of Population Spike (PS) Amplitudes:")
print(f"PS in Alcohol-Preferring Rats: {ps_amplitude_alcohol_preferring_mV} mV")
print(f"PS in Sucrose-Preferring Rats: {ps_amplitude_sucrose_preferring_mV} mV")
print("\nComparison of Slc6a11 Knockdown vs. Baseline:")
print(f"PS after Slc6a11 knockdown in sucrose-preferring rats: {ps_amplitude_shRNA_knockdown_mV} mV")
print(f"PS in alcohol-preferring rats (for comparison): {ps_amplitude_alcohol_preferring_mV} mV")
