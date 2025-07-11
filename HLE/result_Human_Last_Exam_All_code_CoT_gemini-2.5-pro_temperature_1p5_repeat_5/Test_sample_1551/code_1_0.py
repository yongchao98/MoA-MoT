import math

# Part 1: Calculate the expected value of X_19
n_19 = 19
# The formula for odd n is (n-1)*(n-3)/2
val_19 = int((n_19 - 1) * (n_19 - 3) / 2)
# We show the equation and the result.
ans_1 = f"({n_19}-1)*({n_19}-3)/2 = {val_19}"
ans_1_final_val = val_19 # as per question

# Part 2: The expected value of X_20
# For even n, the distance between gifts always has an odd parity,
# but the absorbing states (distance 2 or n-2) are even.
# Thus, the game never ends.
ans_2 = "inf"

# Part 3: General formula for odd n > 1
# This is the formula derived in the thinking steps.
ans_3 = "(n-1)*(n-3)/2"

# Part 4: Expected number of times for distance 10 or n-12 (for n > 30)
# This corresponds to a clockwise distance D of 11 or n-11.
# The formula for the expected number of visits is 2*(n-13).
ans_4 = "2*(n-13)"

# Part 5: Does the game end with probability 1 for every odd n > 1?
# Yes, because for odd n, the random walk on the distance is an irreducible
# Markov chain on a finite state space with absorbing states.
ans_5 = "yes"

# Print the final answers, separated by commas
print(f"{ans_1_final_val}, {ans_2}, {ans_3}, {ans_4}, {ans_5}")
