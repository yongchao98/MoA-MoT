import math

# Step 1: Define the baseline parameters from the Standard Model of Cosmology.
# N_eff_SM is the standard effective number of neutrino species.
N_eff_SM = 3.044

# The total radiation energy density (rho_rad) is related to the photon energy density (rho_gamma) by:
# rho_rad = rho_gamma * (1 + (7/8) * (4/11)^(4/3) * N_eff)
# Let's calculate the conversion factor between N_eff and its contribution to the energy density.
energy_density_factor = (7.0/8.0) * (4.0/11.0)**(4.0/3.0)

# Step 2: Introduce the new physics.
# A new particle decays, injecting extra energy *only* into the neutrino sector.
# Let's assume this injected energy density, delta_rho_nu, is 5% of the photon energy density (rho_gamma).
# This is a hypothetical but illustrative value.
delta_rho_nu_over_rho_gamma = 0.05

# Step 3: Calculate the change in N_eff.
# The additional energy density (delta_rho_nu) contributes to N_eff.
# The change in N_eff is this additional energy density scaled by the conversion factor.
# delta_N_eff = delta_rho_nu / (rho_gamma * energy_density_factor)
delta_N_eff = delta_rho_nu_over_rho_gamma / energy_density_factor

# Step 4: Calculate the new N_eff.
# The new N_eff is the standard model value plus the increase from the new particle's decay.
N_eff_new = N_eff_SM + delta_N_eff

# Step 5: Print the results clearly.
print("--- Analyzing the Impact of a New Decaying Particle on N_eff ---")
print(f"Standard Model N_eff (N_eff_SM): {N_eff_SM:.4f}")
print(f"Assumed energy injected by new particle (as a fraction of photon energy): {delta_rho_nu_over_rho_gamma}")
print("\nThe decay of the new particle adds energy to the neutrino background.")
print("This additional energy density increases the total 'dark radiation', which by definition increases N_eff.")
print("\n--- Calculation ---")
print(f"The increase in N_eff is calculated as: {delta_N_eff:.4f}")
print("The new N_eff is the sum of the old value and the increase.")
print("Final Equation with numbers:")
print(f"{N_eff_new:.4f} = {N_eff_SM:.4f} + {delta_N_eff:.4f}")

print("\nConclusion: The value of N_eff would increase.")
