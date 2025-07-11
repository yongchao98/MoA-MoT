import sys

def count_grid_squares(n_str):
    """
    Calculates the total number of squares on an n x n grid by evaluating
    the summation Sum_{m=1 to n} a^2 * b, where a=n-m+1 and b=m.
    It prints each term of the sum and the final result.

    Args:
        n_str: The size 'n' of the grid, as a string.
    """
    try:
        n = int(n_str)
        if n <= 0:
            print("Error: n must be a positive integer.")
            return
    except ValueError:
        print(f"Error: Invalid input '{n_str}'. Please provide a positive integer for n.")
        return

    total_squares = 0
    terms_list = []
    
    # The formula is Sum_{m=1 to n} (n-m+1)^2 * m
    for m in range(1, n + 1):
        a = n - m + 1
        b = m
        term_value = (a**2) * b
        terms_list.append(str(term_value))
        total_squares += term_value

    # The prompt requires printing each number in the final equation.
    # For n=3, the equation is 9 + 8 + 3 = 20. The numbers are 9, 8, 3, and 20.
    # The following line prints this equation.
    final_equation = " + ".join(terms_list) + f" = {total_squares}"
    
    print(f"For an n={n} grid, the expressions are a=n-m+1 and b=m.")
    print(f"The number of squares is the sum of {n} terms:")
    print(final_equation)

if __name__ == '__main__':
    # The script uses a command-line argument for n.
    # A default value is used if no argument is provided.
    if len(sys.argv) == 2:
        count_grid_squares(sys.argv[1])
    else:
        print("No command-line argument for n found. Running with default value n=3.")
        count_grid_squares("3")
        print("\nUsage example: python your_script_name.py 5")
