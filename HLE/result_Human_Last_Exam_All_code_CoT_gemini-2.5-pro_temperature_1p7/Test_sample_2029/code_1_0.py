import math

def combinations(n, k):
    """
    Calculates the number of combinations (n choose k).
    Returns 0 if k < 0 or k > n.
    """
    if k < 0 or k > n:
        return 0
    # math.comb is available in Python 3.8+ and handles large numbers
    return math.comb(n, k)

def solve():
    """
    Calculates the number of different sequences (a_1,...,a_8, b_1,...,b_6) possible.
    """
    num_a_players = 8
    max_a_wins = 6
    num_b_players = 6
    max_b_wins = 8
    
    total_players = num_a_players + num_b_players
    total_games = num_a_players * num_b_players

    # The calculation is based on the principle of inclusion-exclusion
    # to find the number of integer solutions for a bounded equation.
    total_sequences = 0
    
    # Loop for violations from team A (a_i > 6, so we enforce a_i >= 7)
    for i in range(num_a_players + 1):
        # Loop for violations from team B (b_j > 8, so we enforce b_j >= 9)
        for j in range(num_b_players + 1):
            
            # This is the remaining sum after assuming i players from A and j from B
            # have violated their upper win constraints.
            # Violation for a_i > 6 means a_i >= 7.
            # Violation for b_j > 8 means b_j >= 9.
            remaining_sum = total_games - i * (max_a_wins + 1) - j * (max_b_wins + 1)
            
            if remaining_sum < 0:
                continue
            
            # Number of non-negative solutions for the new sum using stars and bars:
            # C(N + M - 1, M - 1) where N is remaining_sum and M is total_players.
            num_solutions = combinations(remaining_sum + total_players - 1, total_players - 1)
            
            # Number of ways to choose which players violate the constraints.
            ways_to_choose_violators = combinations(num_a_players, i) * combinations(num_b_players, j)
            
            # Apply the inclusion-exclusion sign.
            sign = (-1)**(i + j)
            
            term = sign * ways_to_choose_violators * num_solutions
            total_sequences += term

    # ---- Output the explanation and result ----
    a_vars = " + ".join([f"a_{i}" for i in range(1, num_a_players + 1)])
    b_vars = " + ".join([f"b_{j}" for j in range(1, num_b_players + 1)])
    
    print("This problem is equivalent to finding the number of integer solutions to the equation:")
    print(f"{a_vars} + {b_vars} = {total_games}")
    print("\nSubject to the constraints:")
    print(f"0 <= a_i <= {max_a_wins} for i=1,...,{num_a_players}")
    print(f"0 <= b_j <= {max_b_wins} for j=1,...,{num_b_players}")
    
    print(f"\nThe total number of different sequences is:")
    print(total_sequences)

solve()