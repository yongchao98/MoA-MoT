import math

def display_final_formula(n):
    """
    This function displays the closed-form formula for P(n) for a specific integer n.

    The derivation shows that the general closed-form formula is:
    P(n) = (2 * pi)**(n * (n - 1) / 4) / sqrt(n!)
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: The input 'n' must be a positive integer.")
        return

    # STEP 1: Define the components of the formula based on the input n.
    exponent_numerator = n * (n - 1)
    exponent_denominator = 4
    factorial_value = math.factorial(n)

    # STEP 2: Simplify the exponent for better readability.
    # We find the greatest common divisor to reduce the fraction.
    common_divisor = math.gcd(exponent_numerator, exponent_denominator)
    simplified_exp_num = exponent_numerator // common_divisor
    simplified_exp_den = exponent_denominator // common_divisor

    if simplified_exp_den == 1:
        exponent_str = str(simplified_exp_num)
    else:
        exponent_str = f"{simplified_exp_num}/{simplified_exp_den}"

    # STEP 3: Construct and print the final equation string for the given n.
    # This fulfills the requirement of showing the final equation.
    equation_str = f"P({n}) = (2*pi)^({exponent_str}) / sqrt({factorial_value})"
    print(f"The derived closed-form formula for P(n) evaluated at n = {n} is:")
    print(equation_str)
    
    # STEP 4: Output each number in the final equation as requested.
    print("\nThe numerical components of the equation are:")
    print(f"Base of the power term: 2")
    # pi is a constant, represented symbolically.
    print(f"Another factor in the base: pi")
    if "/" in exponent_str:
        num, den = exponent_str.split('/')
        print(f"Numerator of the exponent: {num}")
        print(f"Denominator of the exponent: {den}")
    else:
        print(f"Exponent: {exponent_str}")
    # The value under the square root comes from n!
    print(f"Value under the square root: {factorial_value}")


# You can change the value of n to see the formula for that specific integer.
n_value = 5
display_final_formula(n_value)
