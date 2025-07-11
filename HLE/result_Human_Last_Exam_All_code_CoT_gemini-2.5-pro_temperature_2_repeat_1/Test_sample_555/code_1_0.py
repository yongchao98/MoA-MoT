import math

def calculate_expected_rolls():
    """
    Calculates the minimal expected number of rolls based on the derived formula.
    The formula for the expected number of rolls E is:
    E = (752 / 343) * (2401 / 2400)
    """
    # Numbers from the final equation
    num1 = 752
    den1 = 343
    num2 = 2401
    den2 = 2400

    print("The final equation for the expected number of rolls, E, is derived as:")
    print(f"E = ({num1}/{den1}) * ({num2}/{den2})")

    # Multiply the fractions
    final_numerator = num1 * num2
    final_denominator = den1 * den2
    
    print(f"Multiplying the fractions gives:")
    print(f"E = {final_numerator}/{final_denominator}")

    # Simplify the resulting fraction
    common_divisor = math.gcd(final_numerator, final_denominator)
    
    simplified_numerator = final_numerator // common_divisor
    simplified_denominator = final_denominator // common_divisor

    print("\nSimplifying the fraction, the minimal expected value of rolls is:")
    print(f"E = {simplified_numerator}/{simplified_denominator}")

calculate_expected_rolls()