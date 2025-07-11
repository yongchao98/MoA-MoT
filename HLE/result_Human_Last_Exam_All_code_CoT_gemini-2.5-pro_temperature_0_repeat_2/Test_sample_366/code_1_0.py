import math

def solve_power_fraction():
    """
    Calculates the fraction of power produced from Pu-239 after 35% uranium burnup.
    """
    # Step 1: Define initial conditions and constants from the problem statement.
    # Initial relative number of atoms based on 30% enrichment.
    N_U235_0 = 30.0
    N_U238_0 = 70.0
    N_total_0 = N_U235_0 + N_U238_0

    # Cross-sections in barns (1 barn = 1e-24 cm^2).
    # We can use the values directly as their units will cancel out.
    sigma_a_U235 = 591.0
    sigma_f_U235 = 505.0
    sigma_a_U238 = 2.42
    sigma_a_Pu239 = 973.0
    sigma_f_Pu239 = 698.0

    # Burnup condition: 35% of initial uranium is consumed.
    burnup_fraction = 0.35
    # The target number of remaining uranium atoms.
    N_U_final_target = N_total_0 * (1.0 - burnup_fraction)

    # Step 2: Find the neutron fluence (tau) that results in the specified burnup.
    # We need to solve the equation:
    # N_U235_0*exp(-sigma_a_U235*tau) + N_U238_0*exp(-sigma_a_U238*tau) = N_U_final_target
    # Let's define a function f(tau) = 0 for the root finding.
    def f(tau):
        return (N_U235_0 * math.exp(-sigma_a_U235 * tau) +
                N_U238_0 * math.exp(-sigma_a_U238 * tau) -
                N_U_final_target)

    # Use the bisection method to find the root tau.
    # Establish a search range [a, b] where f(a) and f(b) have opposite signs.
    a, b = 0.0, 0.1
    if f(a) * f(b) >= 0:
        print("Error: Bisection method requires a valid starting bracket.")
        return

    # Iterate to find the root with high precision.
    for _ in range(100): # 100 iterations is sufficient for double precision.
        midpoint = (a + b) / 2
        if f(midpoint) == 0:
            break
        elif f(a) * f(midpoint) < 0:
            b = midpoint
        else:
            a = midpoint
    tau = (a + b) / 2

    # Step 3: Calculate the final number of atoms for each nuclide using the calculated fluence tau.
    # Final number of U-235 atoms
    N_U235_final = N_U235_0 * math.exp(-sigma_a_U235 * tau)

    # Final number of Pu-239 atoms, using the solution to the Bateman equation.
    # N_Pu(t) = (N_U238_0 * sigma_a_U238 / (sigma_a_Pu239 - sigma_a_U238)) *
    #           (exp(-sigma_a_U238*tau) - exp(-sigma_a_Pu239*tau))
    term1 = N_U238_0 * sigma_a_U238 / (sigma_a_Pu239 - sigma_a_U238)
    term2 = math.exp(-sigma_a_U238 * tau) - math.exp(-sigma_a_Pu239 * tau)
    N_Pu239_final = term1 * term2

    # Step 4: Calculate the power contribution from each fissile isotope.
    # Power is proportional to N * sigma_f.
    power_U235 = N_U235_final * sigma_f_U235
    power_Pu239 = N_Pu239_final * sigma_f_Pu239

    # Step 5: Calculate the fraction of power from Pu-239.
    total_power = power_U235 + power_Pu239
    if total_power == 0:
        fraction_from_Pu = 0.0
    else:
        fraction_from_Pu = power_Pu239 / total_power

    # Step 6: Print the results, including the final equation.
    print("After 35% burnup, the relative number of atoms are:")
    print(f"  N_U235 = {N_U235_final:.6e}")
    print(f"  N_Pu239 = {N_Pu239_final:.6f}\n")

    print("The fraction of power from Plutonium-239 is calculated as:")
    print("Fraction = (N_Pu239 * sigma_f_Pu239) / (N_U235 * sigma_f_U235 + N_Pu239 * sigma_f_Pu239)\n")
    print("Plugging in the values:")
    print(f"Fraction = ({N_Pu239_final:.6f} * {sigma_f_Pu239}) / (({N_U235_final:.6e} * {sigma_f_U235}) + ({N_Pu239_final:.6f} * {sigma_f_Pu239}))")
    print(f"Fraction = {power_Pu239:.6f} / ({power_U235:.6f} + {power_Pu239:.6f})")
    print(f"Fraction = {power_Pu239:.6f} / {total_power:.6f}\n")
    print(f"Final Answer: {fraction_from_Pu:.6f}")
    
    # Return the final answer in the required format
    return fraction_from_Pu

# Execute the function and store the result
final_answer = solve_power_fraction()
# The final answer is printed within the function, but we also output it here in the required format.
print(f"\n<<<{final_answer:.4f}>>>")
