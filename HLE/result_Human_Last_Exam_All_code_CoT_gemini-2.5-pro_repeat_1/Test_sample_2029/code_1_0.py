import math

def combinations(n, k):
    """Calculates 'n choose k', returning 0 for invalid inputs."""
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)

def solve():
    """
    Calculates the number of possible score sequences using generating functions.
    """
    total_sequences = 0
    target_exponent = 48
    team_a_players = 8
    team_b_players = 6
    # The exponent in (1-x)^-k is k = team_a_players + team_b_players = 14
    k_choose = team_a_players + team_b_players - 1 # This is k-1 for C(n+k-1, k-1)

    equation_parts = []
    
    # Iterate through the expansion of (1-x^7)^8, indexed by i
    for i in range(team_a_players + 1):
        # Iterate through the expansion of (1-x^9)^6, indexed by j
        for j in range(team_b_players + 1):
            
            # The power of x from the first two terms is 7*i + 9*j.
            # The remaining power must come from the (1-x)^-14 term.
            k = target_exponent - (7 * i + 9 * j)
            
            if k >= 0:
                # Calculate the coefficient for this combination of i, j, k
                coeff1 = combinations(team_a_players, i) # from (1-x^7)^8
                coeff2 = combinations(team_b_players, j) # from (1-x^9)^6
                coeff3 = combinations(k + k_choose, k_choose) # from (1-x)^-14
                
                term = coeff1 * coeff2 * coeff3
                
                if term == 0:
                    continue
                
                # The sign is determined by (-1)^(i+j)
                if (i + j) % 2 == 1:
                    total_sequences -= term
                    equation_parts.append(f"- {coeff1} * {coeff2} * {coeff3}")
                else:
                    total_sequences += term
                    equation_parts.append(f"+ {coeff1} * {coeff2} * {coeff3}")

    print("The final result is obtained from the following calculation:")
    # Clean up the leading '+' sign for the first term
    if equation_parts[0].startswith('+ '):
        equation_parts[0] = equation_parts[0][2:]
    print(' '.join(equation_parts))
    
    print(f"\nThe total number of different sequences is:")
    print(total_sequences)

solve()
<<<1158922651456>>>