import sys

def calculate_c_k(n, k):
    """
    Calculates the Shapley value for player p_k in the given game.

    Args:
        n (int): The total number of players.
        k (int): The index of the player (from 1 to n).
    """
    if not isinstance(n, int) or not isinstance(k, int) or n <= 1 or k < 1 or k > n:
        print("Invalid input. 'n' must be an integer > 1 and 'k' must be an integer between 1 and n.")
        return

    # S1 is the sum of the first n integers
    S1 = n * (n + 1) // 2
    
    # S2 is the sum of the first n squares
    S2 = n * (n + 1) * (2 * n + 1) // 6

    # The formula for c_k is k*S1^3 + k*S1*S2 - k^2*S1^2
    term1 = k * (S1**3)
    term2 = k * S1 * S2
    term3 = -(k**2) * (S1**2)
    
    c_k = term1 + term2 + term3

    print(f"For n = {n} and k = {k}:")
    print(f"Sum of indices S1 = {S1}")
    print(f"Sum of squared indices S2 = {S2}")
    print("\nThe formula for the fair share c_k is: k*S1^3 + k*S1*S2 - k^2*S1^2")
    print("\nPlugging in the numbers:")
    # Using f-string formatting to show the equation
    print(f"c_{k} = {k} * ({S1})^3 + {k} * {S1} * {S2} - ({k})^2 * ({S1})^2")
    print(f"c_{k} = {term1} + {term2} + {term3}")
    print(f"c_{k} = {c_k}")


# --- Main execution ---
# You can change these values to see the result for different n and k
# As an example, let's use n=4 and k=2
try:
    # Example: n=4, k=2. Change as needed.
    n_val = 4
    k_val = 2
    calculate_c_k(n_val, k_val)
except (ValueError, IndexError):
    print("Please provide integer values for n and k. Example: python your_script.py 4 2")
