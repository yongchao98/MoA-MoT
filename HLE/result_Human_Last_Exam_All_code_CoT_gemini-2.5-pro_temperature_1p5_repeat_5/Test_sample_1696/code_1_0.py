import math

def solve_for_p_n(n):
    """
    This function calculates P(n) using the derived closed-form formula and
    prints the formula and its components.
    
    The closed-form formula for P(n) is:
    P(n) = (2 * pi)^(n*(n-1)/4) / sqrt(n!)
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: Input must be a positive integer.")
        return

    # Calculate components of the formula
    two_pi = 2 * math.pi
    exponent_numerator = n * (n - 1)
    exponent = exponent_numerator / 4
    n_factorial = math.factorial(n)
    sqrt_n_factorial = math.sqrt(n_factorial)
    
    # Calculate the final value, using logs for numerical stability
    try:
        log_p_n = exponent * math.log(two_pi) - 0.5 * math.log(n_factorial)
        final_value = math.exp(log_p_n)
    except (ValueError, OverflowError) as e:
        final_value = float('inf') # Handle large results
        print(f"Could not compute the exact value due to numerical limits: {e}")

    # --- Output ---
    print("The derived closed-form formula is:")
    print("P(n) = (2 * pi)^(n*(n-1)/4) / sqrt(n!)")
    print("\nFor n =", n)
    
    print("The components of the equation are:")
    # Using '2' and 'pi' as requested for the expression
    print("Base of the power:", "2 * pi", f"~ {two_pi:.4f}")
    # Breaking down the exponent calculation
    print("Exponent:", f"n*(n-1)/4 = {n}*({n}-1)/4 = {exponent_numerator}/4 =", exponent)
    # Breaking down the denominator calculation
    print("Denominator:", f"sqrt(n!) = sqrt({n_factorial}) ~ {sqrt_n_factorial:.4e}")
    
    print("\nFinal calculated value:")
    print(f"P({n}) =", final_value)

# Example usage with n=10
# You can change this value to calculate P(n) for any positive integer.
n_value = 10
solve_for_p_n(n_value)