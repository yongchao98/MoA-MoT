import math

def solve_parrot_problem():
    """
    This script finds the best way for the parrot to calculate the mass
    and prints the instructions.
    """
    
    # 1. Define the given values as fractions with integers <= 10.
    density_num = 9
    density_den = 10
    
    radius_num = 1
    radius_den = 2
    
    # 2. Define the constants from the volume formula.
    volume_const_num = 4
    volume_const_den = 3
    
    # 3. Find the best approximation for pi.
    # The acceptable range is [2.827, 3.456].
    # Possible fractions (a/b with a,b <= 10): 3/1, 10/3.
    # The parrot prefers the smallest integers, so we choose 3/1.
    pi_approx_num = 3
    pi_approx_den = 1
    
    # 4. Assemble the instruction for the parrot.
    # The parrot needs to see the full calculation with each number.
    # We expand the radius^3 term to be explicit.
    print("To estimate the mass of the rock, please instruct the parrot to perform the following calculation:")
    print("\nMass = (Density) * (Volume Constant) * (pi approximation) * (radius)^3\n")
    print(f"Mass = ({density_num} / {density_den}) * ({volume_const_num} / {volume_const_den}) * ({pi_approx_num} / {pi_approx_den}) * ({radius_num} / {radius_den}) * ({radius_num} / {radius_den}) * ({radius_num} / {radius_den})")
    
    # 5. Calculate the result of this estimation.
    # Multiply all numerators and denominators.
    final_numerator = density_num * volume_const_num * pi_approx_num * (radius_num ** 3)
    final_denominator = density_den * volume_const_den * pi_approx_den * (radius_den ** 3)
    
    # Simplify the resulting fraction using the greatest common divisor (GCD).
    common_divisor = math.gcd(final_numerator, final_denominator)
    simplified_num = final_numerator // common_divisor
    simplified_den = final_denominator // common_divisor
    
    print(f"\nThe simplified result of this calculation is: {simplified_num} / {simplified_den} kg")

solve_parrot_problem()
