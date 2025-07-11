import sys

# It seems that the user wants to get the answers for the question.
# Based on my analysis, I will calculate the values and then format them as requested.

# Question 1: E[X_19]
n_19 = 19
# For odd n, E[X_n] = (n^2 - 1) / 2
ex_19 = (n_19**2 - 1) / 2

# Question 2: E[X_20]
# For even n, the game never ends, so the expected time is infinite.
ex_20 = float('inf')

# Question 3: E[X_n] for odd n > 1
# The formula is (n^2 - 1) / 2.
formula_odd_n = "(n^2-1)/2"

# Question 4: Expected number of times for distance 11 or n-11 for odd n > 30
# The distance d is 11 (10 friends between) or n-11 (n-12 friends between).
# My derived formula for the sum of expected visits is 2n - 24.
formula_visits = "2*n-24" # This is a symbolic answer as requested for general n > 30.

# Question 5: Is it true that by probability of one the game ends in finite amount of time for every odd n>1 value?
# My analysis says yes.
prob_ends = "Yes"

# Now, let's assemble the final output string.
# We need to print the values for E[X_19], E[X_20], the formula for E[X_n],
# the formula for the visits, and the Yes/No answer.
# The user wants the final equation for odd n, and for the visits.
# I'll represent infinity with 'inf' as is common.
# I need to output the formula for E[X_n] showing the calculation for n.
# And also for the visits.

n_str = 'n'
e_xn_str = f"({n_str}^2-1)/2"

# For the visits part, the question is "what is the expected number of times".
# It can be interpreted as a symbolic answer in terms of n.
# The number of friends between gifts are 10 (d=11) and n-12 (d=n-11).
# Expected visits to d=11 is n-11.
# Expected visits to d=n-11 is n-13.
# Total is (n-11) + (n-13) = 2n - 24.
visits_str = f"(n-11) + (n-13) = 2*n-24"

# The final output needs to be comma-separated.
# E[X_19], E[X_20], E[X_n] for odd n, visits for odd n>30, prob_ends
# The question "Give these answers without any time units"
# "Remember in the final code you still need to output each number in the final equation!"
# This is a bit ambiguous. I'll interpret it as I should show the formula.

# Outputting the answer in a single line as requested by the user prompt
# E[X_19]
n = 19
ex_19_val = (n**2 - 1) // 2
ex_19_str = f"({n}^2-1)/2 = {ex_19_val}"

# E[X_20]
ex_20_str = "inf"

# E[X_n] for odd n
# The prompt asks me to give the answers, so I should provide the formulas as part of the answer string.
e_xn_formula_str = "(n^2-1)/2"

# Visits for odd n > 30
visits_formula_str = "2*n-24"

final_answer_string = f"{ex_19_val}, {ex_20_str}, {e_xn_formula_str}, {visits_formula_str}, {prob_ends}"

print(final_answer_string)
