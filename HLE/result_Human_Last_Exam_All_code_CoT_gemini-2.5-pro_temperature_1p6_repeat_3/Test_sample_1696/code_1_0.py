import math

def calculate_p_n(n):
    """
    Calculates the value of P(n) using the derived closed-form formula.

    The closed-form formula for P(n) is (2*pi)^(n*(n-1)/4) / sqrt(n!).
    This function prints the steps of the calculation for a given n and the final result.
    """
    if not isinstance(n, int) or n <= 0:
        print("Please provide a positive integer for n.")
        return

    # Part 1: State the derived closed-form formula
    print("Through simplification, the closed-form formula for P(n) is found to be:")
    print("P(n) = (2 * pi)^(n * (n - 1) / 4) / sqrt(n!)")
    print("-" * 30)

    # Part 2: Calculate the components of the formula for the given n
    exponent = n * (n - 1) / 4.0
    
    # Use logs for factorial to handle large n without intermediate overflow
    log_n_factorial = math.lgamma(n + 1) # log(n!)
    
    # We will also try to get the integer value for printing, if possible
    n_factorial_str = f"exp({log_n_factorial:.4f})"
    if n < 171: # math.factorial overflows for n >= 171 on many systems
        n_factorial_val = math.factorial(n)
        n_factorial_str = str(n_factorial_val)

    # Part 3: Print the equation with the specific numbers plugged in
    print(f"For n = {n}, we calculate the values in the formula:")
    print(f"The final equation is P({n}) = (2 * pi)^({n}*({n-1})/4) / sqrt({n}!)")
    print(f"Which evaluates to: P({n}) = (2 * {math.pi:.6f})^({exponent}) / sqrt({n_factorial_str})")
    print("-" * 30)

    # Part 4: Calculate the final numerical value using logarithms for stability
    log_numerator = exponent * math.log(2 * math.pi)
    log_denominator = 0.5 * log_n_factorial
    
    log_p_n = log_numerator - log_denominator
    
    # Part 5: Print the final computed value
    print("The final numerical value is:")
    try:
        result = math.exp(log_p_n)
        print(f"P({n}) = {result}")
    except OverflowError:
        print(f"The result P({n}) is too large to represent as a standard float.")
        print(f"The natural logarithm of the result is: ln(P({n})) = {log_p_n}")


if __name__ == '__main__':
    # Set the value of n for the calculation.
    # The user can change this value.
    n_value = 10
    calculate_p_n(n_value)