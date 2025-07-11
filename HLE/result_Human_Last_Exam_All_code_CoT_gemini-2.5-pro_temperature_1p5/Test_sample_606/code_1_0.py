# In the Standard Model of Cosmology, N_eff has a baseline value.
N_eff_SM = 3.044

# This value corresponds to the energy density of the three known neutrino species.
# For simplicity, we can think of the total relativistic energy density (excluding photons)
# as being proportional to N_eff.
# Let's represent this as:
energy_density_neutrinos_SM = N_eff_SM

print(f"Step 1: The Standard Model baseline for the effective number of neutrino species is N_eff_SM = {N_eff_SM}.")

# A new, heavy particle 'X' exists in the early universe. It's non-relativistic, so it's
# part of the matter density, not the radiation density (and thus doesn't contribute to N_eff initially).
# This particle has a non-negligible abundance and decays, releasing its mass-energy.
# The problem states it decays solely into neutrinos.

# This injected energy adds to the existing neutrino energy density.
# Let's represent this injected energy with a positive value, delta_N_eff.
# The exact value depends on the particle's mass and abundance, but we know it's positive.
# For this example, let's assume it's equivalent to adding 0.5 to N_eff.
delta_N_eff = 0.5

print(f"Step 2: A new particle decays, injecting energy *only* into the neutrino population.")
print(f"         This adds an effective energy equivalent of ΔN_eff = {delta_N_eff}.")

# The new total N_eff will be the sum of the Standard Model contribution and the new contribution
# from the decaying particle.
N_eff_new = N_eff_SM + delta_N_eff

print("\nStep 3: The final N_eff is the sum of the original value and the new injected energy.")
print("The calculation is:")
print(f"N_eff_new = N_eff_SM + ΔN_eff")
print(f"N_eff_new = {N_eff_SM} + {delta_N_eff}")
print(f"N_eff_new = {N_eff_new}")


if N_eff_new > N_eff_SM:
    result = "increase"
else:
    result = "decrease"

print(f"\nConclusion: Since the new value ({N_eff_new}) is greater than the standard model value ({N_eff_SM}), N_eff would {result}.")