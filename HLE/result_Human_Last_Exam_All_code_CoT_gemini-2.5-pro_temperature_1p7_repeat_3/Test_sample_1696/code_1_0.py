import math

def solve_p_n(n):
    """
    Calculates the value of P(n) for a positive integer n using its closed-form formula.
    
    The closed-form formula for P(n) is:
    P(n) = (2*pi)^(n*(n-1)/4) / sqrt(n!)
    
    This function prints the formula with the specific numbers for the given n and the final result.
    """
    if not isinstance(n, int) or n <= 0:
        print("Please provide a positive integer for n.")
        return

    # Calculate the components of the formula
    exponent = n * (n - 1) / 4.0
    base = 2 * math.pi
    
    try:
        n_factorial = math.factorial(n)
        sqrt_n_factorial = math.sqrt(n_factorial)
    except OverflowError:
        print(f"n! for n={n} is too large to compute directly.")
        # We can still compute the result using logarithms
        n_factorial_str = "large number"
        sqrt_n_factorial_str = "large number"
    else:
        n_factorial_str = str(n_factorial)
        sqrt_n_factorial_str = f"{sqrt_n_factorial:.6f}"
        
    # Print the formula with the calculated components
    print(f"The closed-form formula for P(n) has been found as: P(n) = (2*pi)^(n*(n-1)/4) / sqrt(n!)")
    print("-" * 20)
    print(f"For n = {n}, the numbers in the final equation are:")
    print(f"Base of the power: 2 * pi = {base:.6f}")
    print(f"Exponent: n*(n-1)/4 = {exponent}")
    print(f"Factorial term: n! = {n_factorial_str}")

    # Calculate the final value.
    # To avoid overflow/underflow with large numbers, we compute using logarithms.
    # log(P(n)) = (n*(n-1)/4) * log(2*pi) - 0.5 * log(n!)
    log_gamma_n = math.lgamma(n + 1)  # lgamma(n+1) computes log(n!)
    log_p_n = exponent * math.log(base) - 0.5 * log_gamma_n
    
    try:
        result = math.exp(log_p_n)
        print("\nFinal Result:")
        print(f"P({n}) = {result}")
    except OverflowError:
        print("\nFinal Result:")
        print(f"The result for P({n}) is too large to represent as a standard floating-point number.")
        print(f"The logarithm of the result is log(P({n})) = {log_p_n}")


# Example for n=5
n = 5
solve_p_n(n)