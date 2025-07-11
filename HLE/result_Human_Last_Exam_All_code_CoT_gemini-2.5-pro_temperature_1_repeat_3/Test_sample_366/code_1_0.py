import math

def solve_power_fraction():
    """
    Calculates the fraction of power from Pu-239 after a specific fuel burnup.
    """
    # 1. Initial conditions and constants
    initial_enrichment = 0.30
    burnup_fraction = 0.35

    # Normalize initial atom counts to 1.0 total Uranium atoms
    N_U235_0 = initial_enrichment
    N_U238_0 = 1.0 - initial_enrichment
    N_U_total_0 = N_U235_0 + N_U238_0
    
    target_burnup_abs = burnup_fraction * N_U_total_0

    # Cross-sections in barns
    s_a_U235 = 591.0   # absorption for U-235
    s_a_U238 = 2.42    # absorption for U-238
    s_a_Pu239 = 973.0  # absorption for Pu-239
    s_f_U235 = 505.0   # fission for U-235
    s_f_Pu239 = 698.0  # fission for Pu-239

    # 2. Define the function to find the root for tau
    # The function represents the difference between the current burnup and the target burnup.
    # We want to find tau such that this function is zero.
    def burnup_equation(tau):
        # Uranium atoms burned
        u235_burned = N_U235_0 * (1 - math.exp(-s_a_U235 * tau))
        u238_burned = N_U238_0 * (1 - math.exp(-s_a_U238 * tau))
        total_burned = u235_burned + u238_burned
        return total_burned - target_burnup_abs

    # 3. Solve for tau using the bisection method
    low_tau = 0.0
    high_tau = 0.1 
    # Check if the root is within the initial bracket
    if burnup_equation(low_tau) * burnup_equation(high_tau) >= 0:
        print("Error: Root is not bracketed. Adjust low_tau and high_tau.")
        return

    for _ in range(100):  # 100 iterations for high precision
        mid_tau = (low_tau + high_tau) / 2
        if burnup_equation(low_tau) * burnup_equation(mid_tau) < 0:
            high_tau = mid_tau
        else:
            low_tau = mid_tau
    
    tau = (low_tau + high_tau) / 2

    # 4. Calculate final atom numbers using the solved tau
    N_U235_final = N_U235_0 * math.exp(-s_a_U235 * tau)
    
    # The analytical solution for N_Pu(tau) considering its production and destruction
    # N_Pu(t) = (lambda_U238_abs * N_U238(t)) / (lambda_Pu239_abs - lambda_U238_abs) * [exp(-lambda_U238_abs*t) - exp(-lambda_Pu239_abs*t)]
    # where lambda = sigma * phi, and tau = phi * t
    # This simplifies to:
    # N_Pu(tau) = (s_a_U238 * N_U238_0 / (s_a_Pu239 - s_a_U238)) * [exp(-s_a_U238*tau) - exp(-s_a_Pu239*tau)]
    
    term1 = (s_a_U238 * N_U238_0) / (s_a_Pu239 - s_a_U238)
    term2 = math.exp(-s_a_U238 * tau) - math.exp(-s_a_Pu239 * tau)
    N_Pu_final = term1 * term2

    # 5. Calculate power from each isotope and the final fraction
    power_U235 = N_U235_final * s_f_U235
    power_Pu239 = N_Pu_final * s_f_Pu239
    
    if (power_U235 + power_Pu239) == 0:
        fraction = 0.0
    else:
        fraction = power_Pu239 / (power_U235 + power_Pu239)

    # Print the final calculation steps
    print("Final state atom numbers (relative):")
    print(f"  N_U235 = {N_U235_final:.6g}")
    print(f"  N_Pu239 = {N_Pu_final:.6g}\n")
    print("Power is proportional to N * sigma_f. The fraction from Plutonium is:")
    print(f"Fraction = (N_Pu239 * sigma_f_Pu239) / (N_U235 * sigma_f_U235 + N_Pu239 * sigma_f_Pu239)")
    print(f"Fraction = ({N_Pu_final:.6g} * {s_f_Pu239}) / (({N_U235_final:.6g} * {s_f_U235}) + ({N_Pu_final:.6g} * {s_f_Pu239}))")
    print(f"Fraction = {power_Pu239:.6g} / ({power_U235:.6g} + {power_Pu239:.6g})")
    print(f"Fraction = {fraction:.4f}")
    
    # Final answer in the required format
    print(f"\n<<<{fraction:.4f}>>>")

solve_power_fraction()