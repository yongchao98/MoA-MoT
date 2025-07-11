import math

def combinations(n, k):
    """
    Calculates the number of combinations 'n choose k'.
    Returns 0 if k < 0 or k > n.
    """
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)

def solve():
    """
    Calculates the number of different sequences (a_1,...,a_8, b_1,...,b_6) that can occur.
    """
    # Total number of games, which is the required sum of all wins.
    total_wins = 48
    
    # Parameters for team A
    num_a_players = 8
    games_per_a_player = 6
    
    # Parameters for team B
    num_b_players = 6
    games_per_b_player = 8
    
    # The calculation is based on finding the coefficient of x^48 in the expansion of
    # P(x) = (1+x+...+x^6)^8 * (1+x+...+x^8)^6
    # This can be written as (1-x^7)^8 * (1-x^9)^6 * (1-x)^(-14)
    # The formula for the coefficient is:
    # Sum_{i,j} [ (-1)^(i+j) * C(8,i) * C(6,j) * C(48 - 7*i - 9*j + 14 - 1, 14 - 1) ]
    # where C(n,k) is n choose k.
    
    total_sequences = 0
    
    # The exponent in (1-x)^(-k) becomes k-1 in the combination C(n+k-1, k-1)
    d = num_a_players + num_b_players # This is 14
    
    print("The total number of sequences is the result of the following sum:")
    print("----------------------------------------------------------------")

    # Iterate over powers i for (1-x^7)^8 and j for (1-x^9)^6
    for j in range(num_b_players + 1):
        for i in range(num_a_players + 1):
            
            exponent = 7 * i + 9 * j
            if exponent <= total_wins:
                
                # Calculate the coefficient from (1-x)^(-14)
                comb_k = total_wins - exponent
                term_comb = combinations(comb_k + d - 1, d - 1)
                
                # Calculate the full term in the sum
                term = ((-1)**(i + j) *
                        combinations(num_a_players, i) *
                        combinations(num_b_players, j) *
                        term_comb)

                if term != 0:
                    sign = "+" if term > 0 else "-"
                    print(f"Term (i={i}, j={j}): {sign} {abs(term)}")

                total_sequences += term

    print("----------------------------------------------------------------")
    print(f"The final number of different sequences is: {total_sequences}")

solve()