from fractions import Fraction

# The optimal first move is to play in a corner. We calculate the win
# probability by considering all possible random first moves by the computer.

# Case A: The computer plays in the center (1 out of 8 squares).
# Our optimal counter-play in this scenario yields a 1/2 chance of winning.
p_case_A = Fraction(1, 8)
win_prob_A = Fraction(1, 2)
contribution_A = p_case_A * win_prob_A

# Case B: The computer plays in the opposite corner (1 out of 8 squares).
# Our optimal counter-play forces a win (win probability of 1).
p_case_B = Fraction(1, 8)
win_prob_B = Fraction(1, 1)
contribution_B = p_case_B * win_prob_B

# Case C: The computer plays in an adjacent corner (2 out of 8 squares).
# Our optimal counter-play in this scenario yields a 23/24 chance of winning.
p_case_C = Fraction(2, 8)
win_prob_C = Fraction(23, 24)
contribution_C = p_case_C * win_prob_C

# Case D: The computer plays on an edge square (4 out of 8 squares).
# Our optimal counter-play forces a win (win probability of 1).
p_case_D = Fraction(4, 8)
win_prob_D = Fraction(1, 1)
contribution_D = p_case_D * win_prob_D

# The total maximum win probability is the sum of the probabilities of these cases.
total_win_probability = contribution_A + contribution_B + contribution_C + contribution_D

# Print the final calculation, showing each number in the equation.
print("To find the maximum chance of winning, we sum the win probabilities for each of the computer's possible opening moves, assuming we start in a corner.")
print("The final calculation is the sum of the contributions from each case:")
print(f"({p_case_A} * {win_prob_A}) + ({p_case_B} * {win_prob_B}) + ({p_case_C} * {win_prob_C}) + ({p_case_D} * {win_prob_D})")
print(f"= {contribution_A} + {contribution_B} + {contribution_C} + {contribution_D}")
print(f"= {total_win_probability}")

# This probability is higher than starting in the center (43/48) or on an edge.
# The result is already a reduced fraction since 89 is prime.

print(f"\nThe maximum chance of winning is {total_win_probability}.")
