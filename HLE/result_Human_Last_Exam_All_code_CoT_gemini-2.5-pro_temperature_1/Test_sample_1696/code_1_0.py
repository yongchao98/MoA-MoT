import math

def solve_for_p(n):
    """
    This function calculates P(n) based on the derived closed-form formula:
    P(n) = (2*pi)^(n*(n-1)/4) / sqrt(n!)
    
    It prints the formula, the components with the given n, and the final result.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.")
        return

    print(f"Finding the value of P(n) for n = {n}")
    print(f"The closed-form formula is: P(n) = (2*pi)^(n*(n-1)/4) / sqrt(n!)")
    print("-" * 30)

    # Calculate the components of the formula
    exponent_numerator = n * (n - 1)
    exponent = exponent_numerator / 4.0
    
    # Handle potentially large n where factorial overflows
    try:
        n_factorial = math.factorial(n)
        sqrt_n_factorial = math.sqrt(n_factorial)

        # Print the equation with numbers substituted
        print("Substituting n into the formula:")
        print(f"The exponent is n*(n-1)/4 = {n}*({n-1})/4 = {exponent_numerator}/4 = {exponent}")
        print(f"The denominator term is n! = {n}! = {n_factorial}")

        print(f"\nFinal Equation: P({n}) = (2*pi)^({exponent}) / sqrt({n_factorial})")

        # Calculate the final value
        value = (2 * math.pi)**exponent / sqrt_n_factorial
        print(f"\nThe calculated value of P({n}) is: {value}")

    except (OverflowError, ValueError):
        # For large n, use logarithms to calculate the value
        print("n is large, calculating with logarithms to avoid overflow.")
        log_2_pi = math.log(2 * math.pi)
        log_gamma_n_plus_1 = math.lgamma(n + 1) # log(n!)
        
        log_p_n = exponent * log_2_pi - 0.5 * log_gamma_n_plus_1
        
        print("Substituting n into the log-formula: log(P(n)) = (n*(n-1)/4)*log(2*pi) - 0.5*log(n!)")
        print(f"The exponent is n*(n-1)/4 = {n}*({n-1})/4 = {exponent}")
        print(f"log(n!) is approximately {log_gamma_n_plus_1}")
        
        print(f"\nlog(P({n})) = {exponent} * {log_2_pi:.4f} - 0.5 * {log_gamma_n_plus_1:.4f} = {log_p_n:.4f}")
        
        # The value might be too large to represent, but we can try
        try:
            value = math.exp(log_p_n)
            print(f"\nThe calculated value of P({n}) is approximately: {value}")
        except OverflowError:
            print(f"\nThe value of P({n}) is too large to represent as a standard float.")

# --- Execute with an example value for n ---
# You can change this value to any positive integer.
n = 5
solve_for_p(n)
