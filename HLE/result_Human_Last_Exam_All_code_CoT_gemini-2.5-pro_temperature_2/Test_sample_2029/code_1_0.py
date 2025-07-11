import math

def calculate_sequences():
    """
    Calculates the number of different score sequences possible.
    """
    
    # Team A has 8 members, Team B has 6 members.
    team_a_size = 8
    team_b_size = 6
    
    # Each member of team A plays 6 games, each member of team B plays 8 games.
    # Total games = 8 * 6 = 48.
    # a_i is win count for player i in team A, 0 <= a_i <= 6.
    # b_j is win count for player j in team B, 0 <= b_j <= 8.
    # Sum(a_i) + Sum(b_j) = 48.
    
    # We are calculating the coefficient of x^48 in the expansion of:
    # (1 + x + ... + x^6)^8 * (1 + x + ... + x^8)^6
    # This equals the coefficient of x^48 in (1-x^7)^8 * (1-x^9)^6 * (1-x)^(-14).
    
    # The formula for the coefficient is:
    # Sum_{i=0 to 8, j=0 to 6} [ C(8,i)*(-1)^i * C(6,j)*(-1)^j * C(48 - 7i - 9j + 14-1, 14-1) ]
    # where C(n, k) is the binomial coefficient "n choose k".

    print("The final calculation is based on the following formula derived from generating functions:")
    # Using the variable names to show how they fit in the equation.
    # The numbers in the equation are:
    # a_size = 8, b_size = 6, total_games = 48
    # max_win_a = 6 -> exponent is 7
    # max_win_b = 8 -> exponent is 9
    # total_players = 14 -> denominator power is 14
    
    print("N = Sum_{i=0..8, j=0..6} [ C(8,i) * (-1)^i * C(6,j) * (-1)^j * C(48 - 7*i - 9*j + 13, 13) ]\n")

    def combinations(n, k):
        if k < 0 or k > n:
            return 0
        return math.comb(n, k)

    total_count = 0
    # i corresponds to the expansion of (1-x^7)^8
    for i in range(team_a_size + 1):
        # j corresponds to the expansion of (1-x^9)^6
        for j in range(team_b_size + 1):
            
            exponent_sum = 7 * i + 9 * j
            if exponent_sum <= 48:
                
                term_a = combinations(team_a_size, i) * ((-1)**i)
                term_b = combinations(team_b_size, j) * ((-1)**j)
                
                # The exponent k for (1-x)^(-14) must be 48 - 7i - 9j
                k = 48 - exponent_sum
                
                # The coefficient of x^k in (1-x)^(-14) is C(k + 14 - 1, 14 - 1)
                term_c = combinations(k + team_a_size + team_b_size - 1, team_a_size + team_b_size - 1)
                
                total_count += term_a * term_b * term_c
                
    print(f"The total number of different sequences is: {total_count}")

calculate_sequences()