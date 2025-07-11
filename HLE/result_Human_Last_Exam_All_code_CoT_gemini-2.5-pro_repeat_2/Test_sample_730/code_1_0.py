import math

def solve_parrot_math():
    """
    This script demonstrates how a parrot could estimate the mass of a rock
    using only fractions with integers 10 or smaller.
    """
    
    # 1. Define the parameters of the problem as fractions
    density_num, density_den = 9, 10
    radius_num, radius_den = 1, 2
    
    # 2. Define constants from the volume formula (V = 4/3 * pi * r^3)
    vol_const_num, vol_const_den = 4, 3
    
    # 3. Approximate pi with a simple integer to meet the parrot's constraints
    pi_approx_num, pi_approx_den = 3, 1
    
    print("Yes, the parrot can calculate the mass. Here is the calculation:")
    print("-----------------------------------------------------------------")
    print("The formula for mass is: Mass = Density * (4/3) * π * Radius^3")
    print("\nWe are given:")
    print(f"Density = {density_num}/{density_den}")
    print(f"Radius = {radius_num}/{radius_den}")
    print("\nWe will use a simple approximation for π:")
    print(f"π ≈ {pi_approx_num}")
    
    print("\nThe full calculation is:")
    print(f"Mass ≈ ({density_num}/{density_den}) * ({vol_const_num}/{vol_const_den}) * ({pi_approx_num}/{pi_approx_den}) * ({radius_num}/{radius_den})^3")
    
    # To show each number, we expand the radius cubed term
    print("\nBreaking it down, we multiply the fractions:")
    print(f"Mass ≈ ({density_num}/{density_den}) * ({vol_const_num}/{vol_const_den}) * ({pi_approx_num}/{pi_approx_den}) * ({radius_num}/{radius_den}) * ({radius_num}/{radius_den}) * ({radius_num}/{radius_den})")
    
    # 4. Calculate the numerator and denominator of the final fraction
    final_numerator = density_num * vol_const_num * pi_approx_num * (radius_num ** 3)
    final_denominator = density_den * vol_const_den * pi_approx_den * (radius_den ** 3)

    print(f"\nMultiplying the numerators: {density_num} * {vol_const_num} * {pi_approx_num} * {radius_num}* {radius_num}* {radius_num} = {final_numerator}")
    print(f"Multiplying the denominators: {density_den} * {vol_const_den} * {pi_approx_den} * {radius_den} * {radius_den} * {radius_den} = {final_denominator}")
    
    # 5. Simplify the resulting fraction
    common_divisor = math.gcd(final_numerator, final_denominator)
    simplified_num = final_numerator // common_divisor
    simplified_den = final_denominator // common_divisor
    
    print(f"\nThe result is {final_numerator}/{final_denominator}, which simplifies to:")
    print(f"Estimated Mass = {simplified_num}/{simplified_den} kg")
    print("-----------------------------------------------------------------")

solve_parrot_math()