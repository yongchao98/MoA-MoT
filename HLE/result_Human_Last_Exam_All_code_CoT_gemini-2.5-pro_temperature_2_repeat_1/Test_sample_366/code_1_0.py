import math

def solve_burnup_fraction():
    """
    Calculates the fraction of power produced from Pu-239 after 35% of the
    initial uranium has been burned.
    """
    # --- 1. Initial Conditions and Constants ---
    # Cross-sections in barns
    sigma_a_U235 = 591.0
    sigma_f_U235 = 505.0
    sigma_a_U238 = 2.42
    sigma_a_Pu239 = 973.0
    sigma_f_Pu239 = 698.0

    # Initial relative number of atoms based on 30% enrichment
    N0_U235 = 30.0
    N0_U238 = 70.0
    N0_total = N0_U235 + N0_U238

    # Burnup condition: 35% of uranium is consumed
    burnup_fraction = 0.35
    N_U_final_target = N0_total * (1.0 - burnup_fraction)

    # --- 2. Solve for Neutron Fluence (Psi) ---
    # We need to solve f(psi) = 0 for psi, where:
    # f(psi) = N0_U235*exp(-sigma_a_U235*psi) + N0_U238*exp(-sigma_a_U238*psi) - N_U_final_target
    
    def f_psi(psi):
        return (N0_U235 * math.exp(-sigma_a_U235 * psi) +
                N0_U238 * math.exp(-sigma_a_U238 * psi) -
                N_U_final_target)

    # Bisection method to find the root psi
    low = 0.0
    high = 0.1 # An initial guess for the upper bound of psi
    tolerance = 1e-12
    
    if f_psi(low) * f_psi(high) >= 0:
        print("Error: Root not bracketed or function does not cross zero.")
        return

    while (high - low) / 2.0 > tolerance:
        mid = (low + high) / 2.0
        if f_psi(mid) == 0:
            psi_final = mid
            break
        elif f_psi(low) * f_psi(mid) < 0:
            high = mid
        else:
            low = mid
    
    psi_final = (low + high) / 2.0

    # --- 3. Calculate Final Atom Counts ---
    # Final number of U-235 atoms
    N_U235_final = N0_U235 * math.exp(-sigma_a_U235 * psi_final)

    # Final number of Pu-239 atoms using Bateman equation solution
    # N_Pu(t) = (lambda_A / (lambda_B - lambda_A)) * N_A(0) * [exp(-lambda_A*t) - exp(-lambda_B*t)]
    # Here, lambda ~ sigma_a, t ~ psi, A=U238, B=Pu239
    term1 = sigma_a_U238 * N0_U238 / (sigma_a_Pu239 - sigma_a_U238)
    term2 = math.exp(-sigma_a_U238 * psi_final) - math.exp(-sigma_a_Pu239 * psi_final)
    N_Pu239_final = term1 * term2

    # --- 4. Calculate Power Fraction ---
    # Power is proportional to N * sigma_f
    power_from_U235 = N_U235_final * sigma_f_U235
    power_from_Pu239 = N_Pu239_final * sigma_f_Pu239
    total_power = power_from_U235 + power_from_Pu239
    
    if total_power == 0:
        fraction_Pu = 0.0
    else:
        fraction_Pu = power_from_Pu239 / total_power

    # --- 5. Print Results ---
    print(f"Solved for total neutron fluence (Ψ): {psi_final:.6g}")
    print(f"Final relative number of U-235 atoms: {N_U235_final:.6g}")
    print(f"Final relative number of Pu-239 atoms: {N_Pu239_final:.6g}\n")
    print("The fraction of power from Pu-239 is calculated by:")
    print("Fraction = (N_Pu239 * σ_f_Pu239) / (N_U235 * σ_f_U235 + N_Pu239 * σ_f_Pu239)\n")
    print("Plugging in the calculated and given values:")
    print(f"Fraction = ({N_Pu239_final:.6g} * {sigma_f_Pu239}) / (({N_U235_final:.6g} * {sigma_f_U235}) + ({N_Pu239_final:.6g} * {sigma_f_Pu239}))")
    print(f"Fraction = {power_from_Pu239:.6g} / ({power_from_U235:.6g} + {power_from_Pu239:.6g})")
    print(f"Fraction = {power_from_Pu239:.6g} / {total_power:.6g}")
    print(f"\nThe fraction of power from Pu-239 is: {fraction_Pu:.5f}")
    
    return fraction_Pu

# Execute the calculation and print the final answer in the desired format
final_answer = solve_burnup_fraction()
# The required output format is <<<value>>>
print(f"\n<<<{final_answer:.5f}>>>")