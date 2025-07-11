import math

def solve_power_fraction():
    """
    Calculates the fraction of power from Pu-239 after 35% of uranium is burned.
    """

    # Step 1: Define initial conditions and cross-sections from the problem
    # Initial atomic fractions (normalized to 1)
    N0_U235 = 0.30
    N0_U238 = 0.70
    N_U_total_initial = N0_U235 + N0_U238

    # Absorption (a) and fission (f) cross-sections in barns
    sigma_a_U235 = 591.0
    sigma_f_U235 = 505.0
    sigma_a_U238 = 2.42
    sigma_a_Pu239 = 973.0
    sigma_f_Pu239 = 698.0
    
    print("Step 1: Initial Conditions and Constants")
    print(f"Initial U-235 fraction (N0_U235): {N0_U235}")
    print(f"Initial U-238 fraction (N0_U238): {N0_U238}")
    print("-" * 30)

    # Step 2: Define the target burnup condition
    # 35% of uranium is burned, so 65% remains.
    final_U_fraction = 1.0 - 0.35
    print("Step 2: Target Burnup Condition")
    print(f"Fraction of Uranium remaining: {final_U_fraction}")
    print("-" * 30)

    # Step 3 & 4: Numerically solve for the required neutron fluence (ψ)
    # We need to find ψ such that: N0_U235*exp(-σa_U235*ψ) + N0_U238*exp(-σa_U238*ψ) = 0.65
    def uranium_remaining_eq(fluence):
        return N0_U235 * math.exp(-sigma_a_U235 * fluence) + \
               N0_U238 * math.exp(-sigma_a_U238 * fluence) - \
               final_U_fraction

    # Use a simple numerical search (bisection method) to find the root
    low_fluence, high_fluence = 0.0, 1.0
    fluence_f = 0.0
    for _ in range(100): # 100 iterations for high precision
        mid_fluence = (low_fluence + high_fluence) / 2
        result = uranium_remaining_eq(mid_fluence)
        if result > 0:
            low_fluence = mid_fluence
        else:
            high_fluence = mid_fluence
    fluence_f = (low_fluence + high_fluence) / 2
    
    print("Step 3 & 4: Solve for Neutron Fluence (ψ)")
    print(f"The fluence (ψ) required for 35% burnup is: {fluence_f:.6f}")
    print("-" * 30)

    # Step 5: Calculate the final concentrations of U-235 and Pu-239 at fluence_f
    # Final concentration of U-235
    Nf_U235 = N0_U235 * math.exp(-sigma_a_U235 * fluence_f)

    # Final concentration of Pu-239 using the Bateman equation for a parent-daughter chain
    # N_Pu(t) = [λ_parent / (λ_daughter - λ_parent)] * N_parent(0) * [exp(-λ_parent*t) - exp(-λ_daughter*t)]
    # Here, λ is rate, which is σ*Φ. So λ*t = σ*ψ.
    bateman_const = (N0_U238 * sigma_a_U238) / (sigma_a_Pu239 - sigma_a_U238)
    Nf_Pu239 = bateman_const * (math.exp(-sigma_a_U238 * fluence_f) - math.exp(-sigma_a_Pu239 * fluence_f))

    print("Step 5: Calculate Final Nuclide Concentrations")
    print(f"Final U-235 concentration (Nf_U235): {Nf_U235:.4e}")
    print(f"Final Pu-239 concentration (Nf_Pu239): {Nf_Pu239:.6f}")
    print("-" * 30)

    # Step 6: Calculate power contributions and the final fraction
    # Power is proportional to the fission rate (N * σ_f)
    power_U235 = Nf_U235 * sigma_f_U235
    power_Pu239 = Nf_Pu239 * sigma_f_Pu239

    total_power = power_U235 + power_Pu239
    
    # Handle the case where total power is zero to avoid division by zero
    if total_power == 0:
        fraction_pu = 0
    else:
        fraction_pu = power_Pu239 / total_power
        
    print("Step 6: Calculate Power Fraction from Plutonium-239")
    print(f"Final relative power from U-235 = {Nf_U235:.4e} * {sigma_f_U235} = {power_U235:.4e}")
    print(f"Final relative power from Pu-239 = {Nf_Pu239:.6f} * {sigma_f_Pu239} = {power_Pu239:.4f}")
    print("\nFraction of Power from Pu-239 = Power(Pu-239) / (Power(U-235) + Power(Pu-239))")
    print(f"Fraction = {power_Pu239:.4f} / ({power_U235:.4e} + {power_Pu239:.4f})")
    print(f"Fraction = {fraction_pu:.6f}")
    print("-" * 30)

    return fraction_pu

# Execute the calculation and print the final answer in the requested format
final_fraction = solve_power_fraction()
# The result is extremely close to 1.0 due to the complete burnout of U-235
# Rounding to a reasonable number of significant figures gives 1.0
final_answer = round(final_fraction, 2)
print(f"The final calculated fraction of power from plutonium-239 is {final_answer}")
<<<1.0>>>