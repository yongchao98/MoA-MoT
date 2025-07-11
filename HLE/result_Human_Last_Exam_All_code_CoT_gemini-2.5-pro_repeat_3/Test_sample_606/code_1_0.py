import math

def calculate_neff(rho_nu, rho_gamma):
    """
    Calculates Neff based on the ratio of neutrino to photon energy densities.

    The formula is derived from: rho_rad = rho_gamma * (1 + (7/8)*(T_nu/T_gamma)^4 * Neff)
    This implies: rho_nu = Neff * (7/8)*(T_nu/T_gamma)^4 * rho_gamma
    Therefore: Neff = (rho_nu / rho_gamma) / ((7/8)*(T_nu/T_gamma)^4)
    """
    # In the Standard Model, after e+- annihilation, T_nu/T_gamma = (4/11)^(1/3)
    t_nu_over_t_gamma_ratio = (4.0 / 11.0)**(1.0 / 3.0)
    
    # This is the constant factor C such that Neff = (rho_nu / rho_gamma) / C
    conversion_factor = (7.0 / 8.0) * (t_nu_over_t_gamma_ratio**4)
    
    neff = (rho_nu / rho_gamma) / conversion_factor
    return neff, conversion_factor

# --- 1. Standard Model Baseline ---
# We start with the standard model prediction for Neff.
N_eff_SM = 3.044

# For simplicity, let's assume the photon energy density (rho_gamma) is 100 arbitrary units.
rho_gamma_val = 100.0

# From this, we can calculate the expected neutrino energy density in the Standard Model.
# Rearranging the formula: rho_nu = Neff * C * rho_gamma
_, C = calculate_neff(1, 1) # Calculate the conversion factor
rho_nu_SM = N_eff_SM * C * rho_gamma_val

print("--- Standard Cosmological Model ---")
print(f"Assuming Photon Energy Density (rho_gamma) = {rho_gamma_val:.1f}")
print(f"The Standard Model Neff is {N_eff_SM}")
print(f"This corresponds to a Neutrino Energy Density (rho_nu_SM) of {rho_nu_SM:.4f}")
print("-" * 40)

# --- 2. New Physics Scenario ---
# A new, out-of-equilibrium particle decays, injecting energy *only* into neutrinos.
# Let's say the added energy density from these decays is 10 arbitrary units.
delta_rho_from_X = 10.0

# The new neutrino energy density is the original one plus the injected energy.
rho_nu_new = rho_nu_SM + delta_rho_from_X

# Now, we calculate the new Neff based on this increased neutrino energy density.
N_eff_new, conversion_factor_new = calculate_neff(rho_nu_new, rho_gamma_val)

print("\n--- New Physics Scenario ---")
print(f"A heavy particle decays, injecting {delta_rho_from_X:.1f} units of energy into neutrinos.")
print(f"Photon Energy Density (rho_gamma) remains {rho_gamma_val:.1f}")
print(f"New Neutrino Energy Density is now {rho_nu_SM:.4f} + {delta_rho_from_X:.1f} = {rho_nu_new:.4f}")

# Final equation for Neff_new
print("\nFinal calculation for the new Neff:")
print(f"Neff_new = ({rho_nu_new:.4f} / {rho_gamma_val:.1f}) / {conversion_factor_new:.4f}")
print(f"Resulting new Neff = {N_eff_new:.4f}")
print("-" * 40)

# --- 3. Conclusion ---
print(f"\nComparing the new value ({N_eff_new:.4f}) to the standard one ({N_eff_SM}), we see it has increased.")
print("Therefore, the decay of a heavy particle into neutrinos would increase Neff.")

<<<increase>>>