import math
from fractions import Fraction

def bernoulli(n):
    """
    Computes the n-th Bernoulli number using the Akiyama-Tanigawa algorithm.
    This ensures exact fractional results.
    """
    if n < 0:
        raise ValueError("Bernoulli number index must be non-negative.")
    A = [Fraction(0)] * (n + 1)
    for m in range(n + 1):
        A[m] = Fraction(1, m + 1)
        for j in range(m, 0, -1):
            A[j-1] = j * (A[j-1] - A[j])
    return A[0]

def double_factorial(n):
    """Computes the double factorial n!!."""
    if not isinstance(n, int) or n < 0:
        raise ValueError("Double factorial is only defined for non-negative integers.")
    if n == 0:
        return 1
    return math.prod(range(n, 0, -2))

def calculate_lambda_integral(g):
    """
    Calculates the integral of the product of lambda classes for a given genus g
    using a specific formula valid for g=3.
    """
    
    # Calculate the components of the formula:
    # |B_{2g}| * (2g - 1)!! / (2g)!
    
    n_bernoulli = 2 * g
    b_2g = bernoulli(n_bernoulli)
    abs_b_2g = abs(b_2g)

    n_double_fact = 2 * g - 1
    df_val = double_factorial(n_double_fact)

    n_fact = 2 * g
    f_val = math.factorial(n_fact)

    # Print the equation with values substituted
    print(f"The integral is calculated using the formula: |B_{2g}| * (2g - 1)!! / (2g)! for g={g}")
    print(f"= |B_{n_bernoulli}| * {n_double_fact}!! / {n_fact}!")
    print(f"= |{b_2g}| * {df_val} / {f_val}")
    print(f"= {abs_b_2g} * {df_val} / {f_val}")
    
    # Perform the final calculation
    numerator = abs_b_2g.numerator * df_val
    denominator = abs_b_2g.denominator * f_val
    result = Fraction(numerator, denominator)
    
    print(f"= {abs_b_2g * df_val} / {f_val}")
    print(f"= {result.numerator}/{result.denominator}")
    
    return result

# The problem concerns the moduli of stable curves of genus 3.
g = 3

# Calculate and print the result
final_result = calculate_lambda_integral(g)

# The result is provided in the format "a/b"
# The question requires printing the final result from the code.
# The following line prints the final answer in the requested format for clarity,
# which will then be put into the final answer block.
# print(f"\nFinal Answer: {final_result.numerator}/{final_result.denominator}")
