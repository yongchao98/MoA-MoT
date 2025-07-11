import math
from itertools import combinations

def get_user_input():
    """Gets and validates user input for n and k."""
    while True:
        try:
            n_str = input("Enter the number of people (n > 1): ")
            n = int(n_str)
            if n <= 1:
                print("Please enter a number n > 1.")
                continue
            
            k_str = input(f"Enter the index of the person k (from 1 to {n}): ")
            k = int(k_str)

            if not (1 <= k <= n):
                print(f"Please enter a number k between 1 and {n}.")
                continue
            
            return n, k
        except ValueError:
            print("Invalid input. Please enter integers.")

def calculate_shapley_value_component_wise(n, k):
    """
    Calculates the Shapley value for player k in the given game.
    The formula used is c_k = k^2*S2 - 2*k^3*S1 + 2*k^2*S1^2 + 4*k*E[A^3]
    """

    # Calculate S1 and S2
    s1 = n * (n + 1) // 2
    s2 = n * (n + 1) * (2 * n + 1) // 6

    # Calculate E[A^3]
    # A is the sum of indices of players in a predecessor set of player k.
    # A random predecessor set is equivalent to choosing a subset of N\{k} of size s
    # with probability 1/n for each s from 0 to n-1.
    
    other_players = list(range(1, n + 1))
    other_players.remove(k)
    
    sum_of_expected_A3 = 0
    
    # Sum over sizes of predecessor sets
    for s_size in range(n):
        # Find all predecessor sets of this size
        subsets_at_s_size = list(combinations(other_players, s_size))
        num_subsets = len(subsets_at_s_size)
        
        # Calculate the average of A^3 for this size
        avg_A3_at_s_size = 0
        if num_subsets > 0:
            sum_A3_at_s_size = 0
            for subset in subsets_at_s_size:
                a = sum(subset)
                sum_A3_at_s_size += a**3
            avg_A3_at_s_size = sum_A3_at_s_size / num_subsets
        
        sum_of_expected_A3 += avg_A3_at_s_size

    # E[A^3] is the average over all possible predecessor sets
    e_a3 = sum_of_expected_A3 / n

    # Calculate each term of the formula
    term1 = k**2 * s2
    term2 = -2 * k**3 * s1
    term3 = 2 * k**2 * s1**2
    term4 = 4 * k * e_a3
    
    c_k = term1 + term2 + term3 + term4

    # Print the formula with the computed values
    print(f"For n = {n} and player p_{k}:")
    print(f"Sum of indices S1 = {s1}")
    print(f"Sum of squared indices S2 = {s2}")
    print(f"Expected value of predecessor sum cubed E[A^3] = {e_a3:.2f}")
    print("\nThe amount of money for player p_k is c_k, calculated as:")
    print(f"c_{k} = k^2*S2 - 2*k^3*S1 + 2*k^2*S1^2 + 4*k*E[A^3]")
    print(f"c_{k} = ({k}^2 * {s2}) - (2 * {k}^3 * {s1}) + (2 * {k}^2 * {s1}^2) + (4 * {k} * {e_a3:.2f})")
    print(f"c_{k} = {term1} - {-term2} + {term3} + {term4:.2f}")
    print(f"c_{k} = {c_k:.2f}")

def main():
    n, k = get_user_input()
    calculate_shapley_value_component_wise(n, k)

if __name__ == "__main__":
    main()
