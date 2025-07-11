import math

def final_calculation():
    """
    This function demonstrates the parrot's calculation for the rock's mass.
    It prints each step of the calculation using fractions with small integers.
    """
    # Given values as fractions
    density_num, density_den = 9, 10  # 0.9 kg/cm^3
    radius_num, radius_den = 1, 2      # 0.5 cm

    # Constants from the volume formula
    vol_const_num, vol_const_den = 4, 3

    # The chosen approximation for pi
    pi_approx_num = 3

    # The radius is cubed in the volume formula
    radius_cubed_num = radius_num ** 3
    radius_cubed_den = radius_den ** 3

    print("The parrot can calculate the mass with the following steps:")
    # Print the initial formula with all the numbers
    print(f"mass = ({density_num}/{density_den}) * ({vol_const_num}/{vol_const_den}) * {pi_approx_num} * ({radius_num}/{radius_den})^3")
    
    # Print the formula after cubing the radius
    print(f"mass = ({density_num}/{density_den}) * ({vol_const_num}/{vol_const_den}) * {pi_approx_num} * ({radius_cubed_num}/{radius_cubed_den})")

    # Combine the numerators and denominators
    final_numerator = density_num * vol_const_num * pi_approx_num * radius_cubed_num
    final_denominator = density_den * vol_const_den * radius_cubed_den

    print(f"mass = ({density_num} * {vol_const_num} * {pi_approx_num} * {radius_cubed_num}) / ({density_den} * {vol_const_den} * {radius_cubed_den})")
    
    # Print the unsimplified result
    print(f"mass = {final_numerator} / {final_denominator}")

    # Simplify the fraction by finding the greatest common divisor (GCD)
    common_divisor = math.gcd(final_numerator, final_denominator)
    simplified_num = final_numerator // common_divisor
    simplified_den = final_denominator // common_divisor

    # Print the final simplified fractional mass
    print(f"The estimated mass is: {simplified_num}/{simplified_den} kg")

final_calculation()