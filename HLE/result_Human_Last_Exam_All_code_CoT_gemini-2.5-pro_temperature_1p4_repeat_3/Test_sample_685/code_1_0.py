import sys

def analyze_nim_game_probability(n, m):
    """
    Determines if the first player in 2D-Generalized NIM has a winning probability
    strictly greater than 50% on a random n x m binary matrix.

    Args:
        n (int): The number of rows in the matrix.
        m (int): The number of columns in the matrix.
    """
    print(f"Analyzing game for an {n} x {m} matrix.")

    # The condition for the first player to have a winning probability > 50%
    # simplifies to checking if the number of P-positions N_P(n, m)
    # is strictly less than 2^(n*m - 1).
    
    # Our analysis shows this is equivalent to the condition n * m > 1.
    print("This is true if and only if n * m > 1.")
    
    product = n * m
    is_winning_probable = product > 1
    
    print("\nEvaluating the final equation:")
    # The final equation is n * m > 1.
    # Outputting each number in the equation as requested.
    print(f"Equation: {n} * {m} > 1")
    print(f"Result:   {product} > 1 is {is_winning_probable}")

    final_value = 1 if is_winning_probable else 0
    print(f"\nThe value of f({n}, {m}) is {final_value}.")


if __name__ == "__main__":
    # You can provide n and m as command-line arguments.
    # Example: python your_script.py 5 4
    try:
        n_arg = int(sys.argv[1])
        m_arg = int(sys.argv[2])
    except (IndexError, ValueError):
        # Default values if no arguments are provided.
        print("Using default values n=5, m=4. Provide them as command-line arguments for other values.")
        n_arg = 5
        m_arg = 4

    analyze_nim_game_probability(n_arg, m_arg)

    # Example for the edge case n=1, m=1
    print("\n--- Testing the edge case n=1, m=1 ---")
    analyze_nim_game_probability(1, 1)
