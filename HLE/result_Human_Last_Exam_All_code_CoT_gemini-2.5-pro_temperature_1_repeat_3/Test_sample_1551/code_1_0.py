import math

# Question 1: What is the expected value of X_19?
n_19 = 19
# For odd n, the formula is (n^2 - 1) / 2
expected_X_19_val = (n_19**2 - 1) / 2
expected_X_19_eq = f"({n_19}^2 - 1) / 2"

# Question 2: What is the expected value of X_20?
# For even n, the game never ends.
expected_X_20_val = math.inf

# Question 3: In general, for odd n>1 what is the expected value of X_n?
expected_X_n_formula = "(n^2 - 1) / 2"

# Question 4: For odd n>30 what is the expected number of times that between the two gifts there are exactly 10 and n-12 friends?
# This corresponds to the distance between gifts being 11.
# Based on random walk analysis on a line (approximation for large n), the value is 4.
expected_times_dist_11 = 4

# Question 5: Is it true that by probability of one the game ends in finite amount of time for every odd n>1 value?
# Yes, the underlying Markov chain is finite and irreducible, so the absorbing state is reached with probability 1.
game_ends_prob_one = "Yes"

# Print the results in the required format
# The problem asks to output the equation for the first question.
print(f"{expected_X_19_eq} = {expected_X_19_val},{expected_X_20_val},{expected_X_n_formula},{expected_times_dist_11},{game_ends_prob_one}")
