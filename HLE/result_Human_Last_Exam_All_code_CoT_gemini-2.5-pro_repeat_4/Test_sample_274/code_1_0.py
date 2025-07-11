import sys

def calculate_squares_in_grid(n):
    """
    Calculates the total number of squares on an n x n grid and
    prints the step-by-step summation.
    An n x n grid corresponds to (n+1)x(n+1) points.
    """
    if not isinstance(n, int) or n < 1:
        print("Error: n must be a positive integer.", file=sys.stderr)
        return

    # The formula is Sum_{m=1 to n} m^2 * (n - m + 1)
    total_squares = 0
    sum_term_list = []
    calc_term_list = []
    value_term_list = []

    for m in range(1, n + 1):
        # For each m from 1 to n, we identify a and b
        a = m
        b = n - m + 1
        term_value = a**2 * b
        total_squares += term_value

        # Store string parts for the final equation printout
        sum_term_list.append(f"{a}^2 * ({n}-{m}+1)")
        calc_term_list.append(f"{a**2} * {b}")
        value_term_list.append(str(term_value))

    # Print the final equation step by step
    print(f"For an {n}x{n} grid, the number of squares is found with a={a_str} and b={b_str} in the expression Sum_{{m=1}}^{n} a^2 * b.")
    print("The final equation is:")
    print(f"Total = {' + '.join(sum_term_list)}")
    print(f"      = {' + '.join(calc_term_list)}")
    print(f"      = {' + '.join(value_term_list)}")
    print(f"      = {total_squares}")

if __name__ == '__main__':
    # Define a and b as strings for the explanation
    a_str = "m"
    b_str = "n-m+1"
    
    try:
        # Take n from command-line argument if provided
        n_input = int(sys.argv[1]) if len(sys.argv) > 1 else 3
        calculate_squares_in_grid(n_input)
    except (ValueError, IndexError):
        # Default to n=3 if no valid argument is given
        print("Usage: python your_script_name.py [n]")
        print("No valid integer 'n' provided. Running with default n=3.\n")
        calculate_squares_in_grid(3)
