import math

def count_sequences():
    """
    Calculates the number of different score sequences possible in a tournament
    between two teams A and B.
    """
    # Problem parameters
    num_players_A = 8
    num_players_B = 6
    max_wins_A = 6  # Each player from A plays 6 games
    max_wins_B = 8  # Each player from B plays 8 games
    total_games = num_players_A * num_players_B

    # Parameters for the inclusion-exclusion formula
    num_vars = num_players_A + num_players_B
    comb_k = num_vars - 1 # k in C(n, k), which is 8+6-1=13
    
    # The number of sequences is the coefficient of z^48 in the expansion of:
    # G(z) = (1+...+z^6)^8 * (1+...+z^8)^6
    # This is calculated using the inclusion-exclusion principle.
    # The formula is:
    # N = sum_{i=0 to 8} sum_{j=0 to 6} (-1)^{i+j} * C(8, i) * C(6, j) * C(48 - 7i - 9j + 13, 13)

    print("The number of sequences is calculated using the formula:")
    print(f"N = sum_{{i=0 to {num_players_A}}} sum_{{j=0 to {num_players_B}}} (-1)^{{i+j}} * C({num_players_A}, i) * C({num_players_B}, j) * C({total_games} - i*({max_wins_A}+1) - j*({max_wins_B}+1) + {comb_k}, {comb_k})")
    print("-" * 30)
    
    total_sequences = 0
    
    # Loop for i from 0 to 8 (for team A)
    for i in range(num_players_A + 1):
        # Loop for j from 0 to 6 (for team B)
        for j in range(num_players_B + 1):
            
            # This is the 'n' in C(n, k) for the 'stars and bars' part of the formula
            comb_n = total_games - i * (max_wins_A + 1) - j * (max_wins_B + 1) + comb_k
            
            # If comb_n is negative, math.comb(comb_n, comb_k) would raise a ValueError.
            # The term is zero in this case, so we skip it.
            if comb_n < 0:
                continue
                
            sign = (-1)**(i + j)
            
            # Calculate the term in the summation
            term = sign * math.comb(num_players_A, i) * math.comb(num_players_B, j) * math.comb(comb_n, comb_k)
            
            total_sequences += term
            
    print(f"The total number of different sequences is: {total_sequences}")

count_sequences()