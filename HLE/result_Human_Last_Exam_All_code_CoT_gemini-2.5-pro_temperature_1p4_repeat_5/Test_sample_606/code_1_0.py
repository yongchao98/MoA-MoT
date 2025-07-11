# Define the standard model and new physics parameters for our conceptual model.

# In the Standard Model (SM), N_eff is approximately 3.044.
N_eff_SM = 3.044

# Let's represent the energy density of a single standard neutrino species as a reference unit.
# The absolute value doesn't matter for this demonstration, so we can set it to 1.
rho_neutrino_single_unit = 1.0

# The total energy density in the neutrino sector in the SM is proportional to N_eff_SM.
rho_total_nu_SM = N_eff_SM * rho_neutrino_single_unit

print(f"Standard Model (SM) calculations:")
print(f"  - N_eff in SM: {N_eff_SM}")
print(f"  - Total neutrino energy density in SM (in reference units): {rho_total_nu_SM:.4f}")
print("-" * 40)

# Now, let's introduce the new particle. It decays and adds energy to the neutrino sector.
# Since the particle has a "non-negligible abundance," this added energy is greater than zero.
# Let's assume it adds an amount of energy equivalent to 0.5 neutrino species for this example.
delta_rho_X_unit = 0.5

print(f"New Physics (NP) Scenario:")
print(f"  - Energy density injected by new particle decay (in reference units): {delta_rho_X_unit}")

# The new total energy density in the neutrino sector is the SM density plus the new contribution.
rho_total_nu_NP = rho_total_nu_SM + delta_rho_X_unit

# The new N_eff is calculated from this new total energy density.
N_eff_NP = rho_total_nu_NP / rho_neutrino_single_unit

print(f"  - New total neutrino energy density (in reference units): {rho_total_nu_NP:.4f}")
print("-" * 40)

# Final Result and Conclusion
print("Final Calculation:")
print(f"The new N_eff is calculated as the sum of the old N_eff and the contribution from the new particle.")
print(f"The equation is: N_eff_new = N_eff_SM + (Delta_Rho_X / Rho_single_nu)")
print(f"Plugging in the numbers: {N_eff_NP:.4f} = {N_eff_SM} + ({delta_rho_X_unit} / {rho_neutrino_single_unit})")

if N_eff_NP > N_eff_SM:
    conclusion = "increase"
else:
    conclusion = "decrease or stay the same"

print(f"\nConclusion: Since the new particle injects energy into the neutrino sector, N_eff will {conclusion}.")