import sys

def calculate_sum(n):
    """
    Calculates the closed form for the sum:
    S_n = sum_{k=0 to n} (2k+1)^5 * C(2k, k) * C(2n-2k, n-k)
    where C(n, k) is the binomial coefficient.

    The closed form is:
    S_n = (n+1)^2 * (63n^3 + 119n^2 + 54n + 8) * 4^(n-3)
    """
    if not isinstance(n, int) or n < 0:
        print("Error: Please provide a non-negative integer for n.")
        return

    # Polynomial part of the expression
    # P(n) = (n+1)^2 * (63n^3 + 119n^2 + 54n + 8)
    n_plus_1_sq = (n + 1)**2
    cubic_poly = (63 * n**3 + 119 * n**2 + 54 * n + 8)
    
    # The term (n+1)^2 * cubic_poly is always divisible by 8 for n>=0
    poly_val = n_plus_1_sq * cubic_poly
    
    # Calculate result as (P(n)/8) * 4^n to handle the 4^(-3) factor
    # This avoids float arithmetic and maintains precision.
    result = (poly_val // 8) * (4**n)

    print(f"For n = {n}:")
    
    # To demonstrate the formula, let's print the components
    print(f"  The formula is: (n+1)^2 * (63*n^3 + 119*n^2 + 54*n + 8) * 4^(n-3)")
    print(f"  Plugging in n={n}:")
    print(f"  ({n}+1)^2 * (63*{n}^3 + 119*{n}^2 + 54*{n} + 8) * 4^({n}-3)")
    print(f"  = {n_plus_1_sq} * ({63 * n**3} + {119 * n**2} + {54 * n} + 8) * 4^({n-3})")
    print(f"  = {n_plus_1_sq} * {cubic_poly} * 4^({n-3})")
    print(f"  = {poly_val} * {4**(n-3) if n>=3 else '1/'+str(4**(3-n))}")
    print(f"  = {result}")


if __name__ == '__main__':
    # You can change the value of n here
    # Example usage:
    # calculate_sum(1)
    # calculate_sum(2)
    # calculate_sum(10)
    
    if len(sys.argv) > 1:
        try:
            n_val = int(sys.argv[1])
            calculate_sum(n_val)
        except ValueError:
            print("Invalid input. Please provide an integer.")
    else:
        # Default example if no command-line argument is given
        print("Usage: python your_script_name.py <n>")
        print("Running with default example n=1:")
        calculate_sum(1)
