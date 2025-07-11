import sys

def calculate_squares_expression(n):
    """
    Calculates the total number of squares on an n x n grid and
    prints the full expression based on the formula:
    Sum_{m=1 to n} m^2 * (n-m+1)
    """
    if not isinstance(n, int) or n < 1:
        print("Error: n must be a positive integer.", file=sys.stderr)
        return

    # a and b are identified as a=m and b=n-m+1
    # The expression is Sum_{m=1 to n} a^2 * b
    
    term_strings = []
    term_calc_strings = []
    term_values = []
    total_sum = 0

    for m in range(1, n + 1):
        # In the formula Sum a^2*b, a=m and b=n-m+1
        a = m
        b = n - m + 1
        
        term_value = (a**2) * b
        total_sum += term_value
        
        # Build strings for each part of the output equation
        term_strings.append(f"{a}^2*({n}-{m}+1)")
        term_calc_strings.append(f"{a**2}*{b}")
        term_values.append(str(term_value))

    # Print the explanation and step-by-step calculation
    print(f"To find the number of squares on a {n}x{n} grid, we use the expression:")
    print(f"  Sum_{{m=1 to {n}}} a^2 * b")
    print(f"Where a = m and b = n-m+1\n")
    print("The calculation unfolds as follows:")
    
    # Final equation step-by-step
    # 1. The formula with m and n
    print(f"= {' + '.join(term_strings)}")
    # 2. First calculation step
    print(f"= {' + '.join(term_calc_strings)}")
    # 3. Second calculation step
    print(f"= {' + '.join(term_values)}")
    # 4. Final result
    print(f"= {total_sum}")
    print(f"\nSo, there are {total_sum} squares on a {n}x{n} grid.")

# --- Main execution ---
# You can change this value to test with other grid sizes
n_grid_size = 4
calculate_squares_expression(n_grid_size)