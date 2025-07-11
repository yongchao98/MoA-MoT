# The standard cosmological model predicts a value for the effective
# number of relativistic species, N_eff.
N_eff_SM = 3.044

# A new hypothetical particle decays, injecting energy into the primordial plasma.
# The problem states this decay is solely into neutrinos. This means the energy
# density of the neutrino radiation field increases. We can model this increase
# as a positive contribution to N_eff, which we'll call delta_N_eff.
# Let's assume an arbitrary positive value for this energy injection.
delta_N_eff = 0.5

# The new value of N_eff is the original Standard Model value plus the
# additional energy contribution from the new particle's decay.
N_eff_new = N_eff_SM + delta_N_eff

# --- Output ---
print("The decay of a new particle into neutrinos adds to the total energy density of relativistic species.")
print("This results in a direct increase in the value of N_eff.")
print("\nLet's demonstrate with a calculation:")
print(f"Standard Model N_eff: {N_eff_SM}")
print(f"Contribution from new particle decay (delta_N_eff): {delta_N_eff}")

# As requested, here is the final equation showing each number.
print("\n--- Final Equation ---")
print(f"N_eff_new = N_eff_SM + delta_N_eff")
print(f"{N_eff_new:.3f} = {N_eff_SM} + {delta_N_eff}")
print("---")

print(f"\nThe new value, N_eff_new = {N_eff_new:.3f}, is greater than the original value of {N_eff_SM}.")
print("Therefore, N_eff would increase.")
