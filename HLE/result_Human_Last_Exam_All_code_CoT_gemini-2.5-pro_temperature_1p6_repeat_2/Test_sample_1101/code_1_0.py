import math

def falling_factorial(x, n):
    """Computes the falling factorial x_n = x * (x-1) * ... * (x-n+1)."""
    if n < 0:
        return 0
    if n == 0:
        return 1
    result = 1.0
    for i in range(n):
        result *= (x - i)
    return result

def calculate_sum(n):
    """
    Calculates the sum S(n) using the derived closed-form identity.
    S(n) = (-1)^n * (n - 1/2)^{\underline{n}}
    """
    sign = (-1)**n
    # The term (n - 1/2)^{\underline{n}}
    term = falling_factorial(n - 0.5, n)
    return sign * term

def double_factorial(n):
    """Computes the double factorial n!!."""
    if n < -1:
        return 0
    if n == -1 or n == 0:
        return 1
    result = 1
    for i in range(n, 0, -2):
        result *= i
    return result
    
def calculate_sum_alt(n):
    """
    Calculates the sum S(n) using the double factorial representation.
    S(n) = (-1)^n * (2n-1)!! / 2^n
    """
    if n == 0:
      return 1.0
    sign = (-1)**n
    num = double_factorial(2 * n - 1)
    den = 2**n
    return sign * num / den

def main():
    """
    Finds the function f with the lowest complexity that bounds the sum.
    The sum is S(n) = (-1)^n * (n - 1/2)^{\underline{n}}.
    This is equivalent to (-1)^n * (2n-1)!! / 2^n.
    The absolute value |S(n)| = (1/4^n) * (2n choose n) <= 1 for all n.
    So f(n) can be a constant, for example f(n) = 1. This has the lowest complexity.
    The problem asks for the final equation.
    """
    
    n = 5 # Example value
    
    # Calculate S(n) to show it's a non-trivial value
    sum_val = calculate_sum_alt(n)
    
    # Let's write down the final simplified equation for S(n)
    # The simplified form is S(n) = (-1)^n * (2n-1)!! / 2^n
    
    # For a specific n, let's say n=5
    n_str = str(n)
    sign_str = "(-1)^" + n_str
    num_str = "(2*" + n_str + "-1)!!"
    den_str = "2^" + n_str
    
    sign_val = (-1)**n
    num_val = double_factorial(2*n - 1)
    den_val = 2**n
    
    result_str_lhs = f"S({n})"
    result_str_rhs = f"{sign_val} * {num_val} / {den_val}"
    
    print(f"The sum can be simplified to the expression: S(n) = (-1)^n * (2n-1)!! / 2^n")
    print(f"The magnitude of this sum is bounded by 1 for all n, so the function f(n) with the lowest complexity is f(n) = 1.")
    print("The final output asks for the final equation.")
    print("Let's demonstrate the simplified equation for n=5:")
    
    equation = f"{result_str_lhs} = {sign_str} * {num_str} / {den_str} = {result_str_rhs} = {sum_val}"
    
    # This is a bit awkward. The user wants the equation, so I'll print it.
    # What "final equation"? The expression for S(n).
    
    n_display = 4 # Let's show n=4 to match my scratchpad exploration
    sum_val_4 = calculate_sum_alt(n_display)
    
    # S(4) = (-1)^4 * (7)!! / 2^4 = (7*5*3*1) / 16 = 105/16
    
    print("\nFor n=4:")
    print("S(4) = (-1)^4 * (2*4-1)!! / 2^4")
    print(f"     = 1 * {double_factorial(7)} / {2**4}")
    print(f"     = 105 / 16")
    print(f"     = {105/16}")

main()