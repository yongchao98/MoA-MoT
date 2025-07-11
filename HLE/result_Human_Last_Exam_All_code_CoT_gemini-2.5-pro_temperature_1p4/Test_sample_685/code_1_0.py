import sys

def analyze_nim_game_and_compute_f(n, m):
    """
    This function determines the value of f(n, m) by checking if the first player's
    winning probability is strictly greater than 50%. It shows the calculation steps.

    Args:
        n (int): The number of rows in the matrix.
        m (int): The number of columns in the matrix.
    """
    if not isinstance(n, int) or not isinstance(m, int) or n <= 0 or m <= 0:
        print("Error: Dimensions n and m must be positive integers.")
        return

    print(f"Analyzing the game for a {n}x{m} matrix.")
    
    # As derived in the explanation, the condition for the first player to have a
    # winning probability strictly greater than 50% simplifies to n * m > 1.
    
    # Step 1: Calculate the product of the dimensions.
    product = n * m
    
    # Step 2: Evaluate the final inequality. The request to "output each number
    # in the final equation" refers to showing the values in this condition.
    print(f"\nThe final condition to check is: n * m > 1")
    print(f"Substituting the given values: {n} * {m} > 1")
    print(f"Result of the product: {product}")
    
    is_winning_prob_high = product > 1
    print(f"Evaluating the inequality: {product} > 1 is {is_winning_prob_high}")
    
    # The function f(n, m) returns 1 if the condition is true, 0 otherwise.
    result = 1 if is_winning_prob_high else 0
    print(f"\nTherefore, the value of f({n}, {m}) is: {result}")

if __name__ == '__main__':
    # This part allows running the script from the command line,
    # e.g., python your_script.py 5 10
    if len(sys.argv) == 3:
        try:
            n_val = int(sys.argv[1])
            m_val = int(sys.argv[2])
            analyze_nim_game_and_compute_f(n_val, m_val)
        except ValueError:
            print("Error: Please provide valid integer values for n and m.")
    else:
        # If no command line arguments are given, run some examples.
        print("Running demonstration with two examples.")
        print("-" * 30)
        analyze_nim_game_and_compute_f(5, 10)
        print("-" * 30)
        analyze_nim_game_and_compute_f(1, 1)
        print("-" * 30)
        print("To run with your own values, use: python your_script.py <n> <m>")
