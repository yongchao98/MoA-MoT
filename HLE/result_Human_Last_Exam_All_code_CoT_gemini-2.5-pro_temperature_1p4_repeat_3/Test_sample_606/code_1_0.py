# In the Standard Model of cosmology, the effective number of neutrino species
# has a precisely predicted value.
N_eff_SM = 3.044

# Let's assume the new hypothetical particle exists and decays.
# Its decay injects energy into the neutrino background. We can model this
# by quantifying the additional energy as a fraction of the original
# standard model neutrino energy density.
# For this example, let's assume the decay adds energy equivalent to 15%
# of the pre-existing neutrino energy.
energy_injection_ratio = 0.15

# The new N_eff is proportional to the new total neutrino energy density.
# New Neutrino Energy = (Standard Neutrino Energy) + (Injected Energy)
# New Neutrino Energy = (Standard Neutrino Energy) * (1 + energy_injection_ratio)
# Since N_eff is proportional to the neutrino energy, the new N_eff is:
# N_eff_new = N_eff_SM * (1 + energy_injection_ratio)

new_N_eff = N_eff_SM * (1 + energy_injection_ratio)

# --- Output the results ---
print("The decay of a massive particle into neutrinos adds energy to the neutrino background.")
print("This increases the total radiation density relative to the photon density.")
print("As a result, the value of N_eff increases.")
print("\n--- Illustrative Calculation ---")
print(f"Standard Model N_eff (N_eff_SM): {N_eff_SM}")
print(f"Assumed fractional energy injected by new particle decay: {energy_injection_ratio}")
print("\nThe new N_eff is calculated as: N_eff_new = N_eff_SM * (1 + energy_injection_ratio)")
print("\nFinal Equation:")
print(f"{new_N_eff:.4f} = {N_eff_SM} * (1 + {energy_injection_ratio})")
print(f"\nThe new N_eff is {new_N_eff:.4f}, which is an increase from the standard value of {N_eff_SM}.")
