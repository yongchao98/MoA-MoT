import math

def parrot_mass_calculation():
    """
    This script calculates the mass of the rock using fractions with small integers
    suitable for a clever parrot.

    The formula is: mass = density * volume
    Volume of a sphere = (4/3) * pi * radius^3
    """

    # Given values as fractions with integers <= 10
    density_num = 9
    density_den = 10

    radius_num = 1
    radius_den = 2

    # Constants from the volume formula
    vol_const_num = 4
    vol_const_den = 3

    # Approximate pi with a simple integer that keeps the error low.
    # pi_approx = 3/1 is a valid choice.
    # True mass is ~0.4712 kg.
    # Estimated mass with pi=3 is (9/10)*(4/3)*3*(1/8) = 108/240 = 0.45 kg.
    # Error is |0.45 - 0.4712| / 0.4712 = 4.5%, which is < 10%.
    pi_approx_num = 3
    pi_approx_den = 1 # Not shown in final output for simplicity, as it's 3/1

    # The parrot's calculation:
    # mass = (9/10) * (4/3) * 3 * (1/2)^3
    # First, calculate the numerator and denominator of the result
    final_numerator = density_num * vol_const_num * pi_approx_num * (radius_num**3)
    final_denominator = density_den * vol_const_den * pi_approx_den * (radius_den**3)

    # To show the result with the smallest possible integers, simplify the fraction
    common_divisor = math.gcd(final_numerator, final_denominator)
    simplified_num = final_numerator // common_divisor
    simplified_den = final_denominator // common_divisor

    # Print the full equation for the parrot
    print("Here is the calculation to estimate the mass:")
    print(f"mass = ({density_num}/{density_den}) * ({vol_const_num}/{vol_const_den}) * {pi_approx_num} * ({radius_num}/{radius_den})**3")
    print(f"mass = {simplified_num}/{simplified_den} kg")

parrot_mass_calculation()