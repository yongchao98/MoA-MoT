import math

def solve_fluence():
    """
    Calculates the fraction of power produced by Plutonium-239 fission
    at a specific burnup point.
    """

    # --- 1. Define Constants and Initial Conditions ---
    
    # Cross-sections in barns (1 barn = 1e-24 cm^2)
    s_a_U235 = 591.0   # Absorption for U-235
    s_f_U235 = 505.0   # Fission for U-235
    s_a_U238 = 2.42    # Absorption for U-238
    s_a_Pu239 = 973.0  # Absorption for Pu-239
    s_f_Pu239 = 698.0  # Fission for Pu-239

    # Initial fuel composition (normalized to 1 atom total)
    initial_enrichment = 0.30
    N0_U235 = initial_enrichment
    N0_U238 = 1.0 - initial_enrichment
    
    # Burnup condition
    uranium_burnup_fraction = 0.35
    N_U_target = (N0_U235 + N0_U238) * (1.0 - uranium_burnup_fraction)

    # --- 2. Solve for the Fluence (t) ---
    
    # We need to find t such that: N0_U235*exp(-s_a_U235*t) + N0_U238*exp(-s_a_U238*t) = N_U_target
    # This is equivalent to finding the root of f(t) = 0
    def f(t):
        return N0_U235 * math.exp(-s_a_U235 * t) + N0_U238 * math.exp(-s_a_U238 * t) - N_U_target

    # Simple bisection solver
    low_t, high_t = 0.0, 0.1
    tolerance = 1e-9
    while (high_t - low_t) > tolerance:
        mid_t = (low_t + high_t) / 2
        if f(mid_t) * f(low_t) > 0:
            low_t = mid_t
        else:
            high_t = mid_t
    fluence_t = (low_t + high_t) / 2

    # --- 3. Calculate Final Number of Atoms ---

    # Number of U-235 atoms at the final state
    Nf_U235 = N0_U235 * math.exp(-s_a_U235 * fluence_t)

    # Number of Pu-239 atoms at the final state
    # This comes from the solution to the differential equation for Pu-239 generation and burnout
    prefactor = (N0_U238 * s_a_U238) / (s_a_Pu239 - s_a_U238)
    exp_term_U238 = math.exp(-s_a_U238 * fluence_t)
    exp_term_Pu239 = math.exp(-s_a_Pu239 * fluence_t)
    Nf_Pu239 = prefactor * (exp_term_U238 - exp_term_Pu239)

    # --- 4. Calculate Power Fraction ---
    
    # Power is proportional to N * sigma_fission
    power_from_U235 = Nf_U235 * s_f_U235
    power_from_Pu239 = Nf_Pu239 * s_f_Pu239
    total_power = power_from_U235 + power_from_Pu239
    
    if total_power == 0:
        fraction_from_Pu239 = 0
    else:
        fraction_from_Pu239 = power_from_Pu239 / total_power
        
    # --- 5. Print the Final Calculation Details ---
    print("Calculation of the power fraction from Plutonium-239:")
    print("Fraction = (N_Pu239 * sigma_f_Pu239) / (N_U235 * sigma_f_U235 + N_Pu239 * sigma_f_Pu239)")
    print(f"Fraction = ({Nf_Pu239:.6e} * {s_f_Pu239}) / (({Nf_U235:.6e} * {s_f_U235}) + ({Nf_Pu239:.6e} * {s_f_Pu239}))")
    print(f"Fraction = ({power_from_Pu239:.6e}) / ({power_from_U235:.6e} + {power_from_Pu239:.6e})")
    print(f"Fraction = {fraction_from_Pu239:.6f}")

    return fraction_from_Pu239

if __name__ == "__main__":
    result = solve_fluence()
    # The final answer is wrapped in <<<>>>
    print(f"\nFinal Answer: {result}")
    print(f"<<<{result}>>>")