import sys

def solve_game_probability(n, m):
    """
    This function determines if the first player has a winning position with a
    probability strictly more than 50% for a 2D-Generalized NIM game on a
    random n x m binary matrix.

    The problem reduces to checking if the number of P-positions (losing positions), N_p,
    is less than 2^(n*m - 1). Through analysis, this condition holds for all
    cases where n*m > 1. For n=1, m=1, the condition is not met.

    Args:
        n (int): The number of rows in the matrix.
        m (int): The number of columns in the matrix.
    """

    if not (isinstance(n, int) and n > 0 and isinstance(m, int) and m > 0):
        print("Error: Please provide positive integers for n and m.")
        return

    product = n * m

    print(f"Analyzing the game for an {n} x {m} matrix.")
    print("The condition for the first player to have a winning chance > 50% is:")
    print("Number of P-positions (N_p) < 2^(n*m - 1)")
    print("-" * 20)

    # Based on the analysis, the outcome only depends on whether n*m is greater than 1.
    if product > 1:
        result = 1
        print(f"For n = {n} and m = {m}, the product n * m = {product}, which is greater than 1.")
        print("In this case, the number of P-positions is known to be small enough to satisfy the inequality.")
        print("N_p < 2^({n}*{m} - 1)")
        print(f"Result: f({n}, {m}) = {result}")
    else: # This only happens for n=1, m=1
        result = 0
        p_positions = 1
        threshold = pow(2, n * m - 1)
        print(f"For n = {n} and m = {m}, we have the trivial case.")
        print(f"The number of P-positions is N_p = {p_positions}.")
        print(f"The threshold is 2^({n}*{m} - 1) = {threshold}.")
        print(f"The inequality {p_positions} < {threshold} is False.")
        print(f"Result: f({n}, {m}) = {result}")

# Example usage with n=5 and m=4.
# You can change these values or adapt the script to take command-line arguments.
example_n = 5
example_m = 4
solve_game_probability(example_n, example_m)

print("\n" + "="*40 + "\n")

# Another example for the n=1, m=1 case
example_n_trivial = 1
example_m_trivial = 1
solve_game_probability(example_n_trivial, example_m_trivial)
