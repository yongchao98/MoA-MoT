from fractions import Fraction

# For a large m, there are 11 possible values of n that result in a solution,
# provided that n < N. The 11 values are n=0..8, n=m, and n=3m.

# The condition n < N is always satisfied for N -> infinity for n=0..8 and n=m.
# This gives 10 solutions.
num_solutions_always = 10

# The 11th solution is for n=3m, which is only counted if 3m < N, i.e., m < N/3.
num_solutions_conditional = 11

# As N -> infinity, the fraction of m values such that m < N/3 approaches 1/3.
prob_conditional_met = Fraction(1, 3)

# The fraction of m values such that m >= N/3 approaches 2/3.
prob_conditional_not_met = Fraction(2, 3)

# The limit is the weighted average of the number of solutions.
limit_value = num_solutions_conditional * prob_conditional_met + num_solutions_always * prob_conditional_not_met

# We print out the calculation step by step, including each number in the final equation.
# First, the numbers of solutions in the two cases for m.
n1 = num_solutions_conditional
n2 = num_solutions_always
# Then, the probabilities (proportions) of these cases.
p1 = prob_conditional_met
p2 = prob_conditional_not_met

print(f"The limit is calculated as a weighted average based on the value of m relative to N.")
print(f"For m < N/3 (proportion {p1}), there are {n1} solutions for n.")
print(f"For m >= N/3 (proportion {p2}), there are {n2} solutions for n.")
print(f"The limit is the average: {n1} * {p1} + {n2} * {p2} = {limit_value.numerator}/{limit_value.denominator}")
