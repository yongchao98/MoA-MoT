import math

# Part 1: E[X_19]
n_19 = 19
# For odd n, the formula is (n^2 - 1) / 2
ans_19 = (n_19**2 - 1) / 2

# Part 2: E[X_20]
# For even n, the expected time is infinite
ans_20 = "inf"

# Part 3: General formula for odd n>1
ans_odd_n_formula_str = "(n^2-1)/2"

# Part 4: Expected number of visits for n > 30
# This value is constant for large n.
ans_visits = 4

# Part 5: Probability of game ending for odd n
ans_prob_one = "Yes"

# Print the final answers, comma-separated
print(f"{int(ans_19)},{ans_20},{ans_odd_n_formula_str},{ans_visits},{ans_prob_one}")
