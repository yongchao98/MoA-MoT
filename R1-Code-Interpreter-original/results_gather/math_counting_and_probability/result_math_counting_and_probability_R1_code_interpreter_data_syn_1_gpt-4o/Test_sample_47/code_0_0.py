from math import comb

# Probability of Lakers winning a single game
p_lakers_win = 1/3

# Probability of Celtics winning a single game
p_celtics_win = 2/3

# Number of ways Lakers can win 3 out of the first 6 games
ways_lakers_win_3 = comb(6, 3)

# Probability of Lakers winning exactly 3 out of the first 6 games
prob_lakers_win_3_of_6 = ways_lakers_win_3 * (p_lakers_win**3) * (p_celtics_win**3)

# Probability of Lakers winning the 7th game
prob_lakers_win_7th = p_lakers_win

# Total probability of Lakers winning in exactly 7 games
prob_lakers_win_in_7 = prob_lakers_win_3_of_6 * prob_lakers_win_7th

print(prob_lakers_win_in_7)