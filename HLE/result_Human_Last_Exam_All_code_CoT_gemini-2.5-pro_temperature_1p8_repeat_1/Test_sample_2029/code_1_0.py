import math

def comb(n, k):
    """
    Computes the combination 'n choose k'. Returns 0 if k < 0 or k > n.
    """
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)

# The number of players in team A and B
team_a_players = 8
team_b_players = 6
# Total games played
total_games = team_a_players * team_b_players

# Maximum score for a player in team A and B
max_a_score = team_b_players
max_b_score = team_a_players

# The target exponent for the generating function coefficient
target_sum = total_games

# The problem is to find the coefficient of x^target_sum in the expansion of P(x):
# P(x) = (1+x+...+x^max_a_score)^team_a_players * (1+x+...+x^max_b_score)^team_b_players
# P(x) = ((1-x^(max_a_score+1))/(1-x))^team_a_players * ((1-x^(max_b_score+1))/(1-x))^team_b_players
# P(x) = (1-x^7)^8 * (1-x^9)^6 * (1-x)^(-14)
#
# The coefficient is given by the sum:
# Sum_{i=0..8, j=0..6} C(8,i)*(-1)^i * C(6,j)*(-1)^j * C(48 - 7*i - 9*j + 14 - 1, 14-1)
# where C(n,k) is "n choose k".

total_sequences = 0
is_first_term = True

print("The total number of sequences is given by the sum:")
# Loop over the expansion of (1-x^7)^8
for i in range(team_a_players + 1):
    # Loop over the expansion of (1-x^9)^6
    for j in range(team_b_players + 1):
        
        k = target_sum - (max_a_score + 1) * i - (max_b_score + 1) * j
        
        if k >= 0:
            # Coefficient from (1-x)^(-14) expansion
            power_of_neg_x = team_a_players + team_b_players
            term_val = (comb(team_a_players, i) *
                        comb(team_b_players, j) *
                        comb(k + power_of_neg_x - 1, power_of_neg_x - 1))
            
            # Add or subtract the term based on the sign
            if (i + j) % 2 != 0:
                total_sequences -= term_val
                # Print the term from the equation
                print(f" - C({team_a_players},{i})*C({team_b_players},{j})*C({k + power_of_neg_x - 1},{power_of_neg_x - 1})")
            else:
                total_sequences += term_val
                # Print the term from the equation, handling the sign for the first line
                if is_first_term:
                    print(f"   C({team_a_players},{i})*C({team_b_players},{j})*C({k + power_of_neg_x - 1},{power_of_neg_x - 1})")
                    is_first_term = False
                else:
                    print(f" + C({team_a_players},{i})*C({team_b_players},{j})*C({k + power_of_neg_x - 1},{power_of_neg_x - 1})")

print(f"\n= {total_sequences}")
