# In cosmology, the total energy density of radiation is the sum of the energy density
# of photons and the energy density of neutrinos.
# Neff, the effective number of neutrino species, parameterizes the neutrino energy density.

# For simplicity, we can think of the total neutrino energy density as:
# rho_nu_total = Neff * (energy_of_one_neutrino_species)

# 1. Standard Model Value
# In the Standard Model of cosmology, Neff is slightly greater than 3.
N_eff_SM = 3.044

# We can express the total neutrino energy density in the Standard Model in some arbitrary units.
# Let's say the energy density of a single neutrino species is 1 unit.
energy_per_neutrino_species = 1.0
rho_nu_SM = N_eff_SM * energy_per_neutrino_species
print(f"Standard Model Neff: {N_eff_SM}")
print(f"Standard Model neutrino energy density (relative units): {rho_nu_SM:.3f}")
print("-" * 30)

# 2. Introduce the New Particle's Contribution
# A new particle with a non-negligible abundance decays, injecting energy into the neutrinos.
# Let's assume its energy density is equivalent to half a neutrino species.
rho_X_decay = 0.5 * energy_per_neutrino_species
print(f"Energy injected from new particle decay (relative units): {rho_X_decay}")

# 3. Calculate the New Neutrino Energy Density
# This injected energy is added to the existing neutrino energy density.
rho_nu_new = rho_nu_SM + rho_X_decay

# 4. Calculate the New Neff
# The new Neff is the new total neutrino energy density divided by the energy of a single species.
N_eff_new = rho_nu_new / energy_per_neutrino_species

print(f"After decay, the new total neutrino energy density is:")
print(f"{rho_nu_SM:.3f} (from SM neutrinos) + {rho_X_decay} (from new particle) = {rho_nu_new:.3f}")

print("\nThe new value of Neff is calculated from this new energy density:")
print(f"N_eff_new = {rho_nu_new:.3f} / {energy_per_neutrino_species} = {N_eff_new:.3f}")

if N_eff_new > N_eff_SM:
    print("\nConclusion: The value of Neff has increased.")
else:
    print("\nConclusion: The value of Neff has not increased.")
