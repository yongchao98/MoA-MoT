from fractions import Fraction

# Step 1: Define the basic probabilities for a single coin flip.
P_H = Fraction(1, 3)
P_T = Fraction(2, 3)

# Step 2: Calculate the probability of the events needed.
# Event A: "All three coins are heads" (HHH).
# This is also the case for "0 tails".
P_all_heads = P_H**3

# Event B: "The number of tails is even". This means 0 tails or 2 tails.
# Probability of 0 tails is P_all_heads.
P_0_tails = P_all_heads

# Probability of 2 tails.
# There are 3 combinations for 2 tails (HTT, THT, TTH).
# The probability of any single one of these combinations is P(H) * P(T)^2.
P_one_outcome_with_2_tails = P_H * (P_T**2)
num_combinations_for_2_tails = 3
P_2_tails = num_combinations_for_2_tails * P_one_outcome_with_2_tails

# The total probability of the number of tails being even is P(0 tails) + P(2 tails).
P_even_tails = P_0_tails + P_2_tails

# Step 3: Calculate the conditional probability P(A|B).
# P(A|B) = P(A and B) / P(B)
# The event "A and B" (all heads and even tails) is the same as event "A" (all heads),
# because 'all heads' means 0 tails, and 0 is an even number.
# So, P(A and B) = P(A) = P_all_heads.
conditional_probability = P_all_heads / P_even_tails

# Step 4: Print the final result in the form of an equation.
print("We want to find P(All Heads | Even Number of Tails).")
print("This is calculated as P(All Heads) / P(Even Number of Tails).")
print("\nP(All Heads) = P(HHH) = (1/3)^3 =", P_all_heads)
print("\nP(Even Number of Tails) = P(0 tails) + P(2 tails)")
print(f"P(Even Number of Tails) = {P_0_tails} + 3 * (1/3 * (2/3)^2) = {P_0_tails} + {P_2_tails} = {P_even_tails}")
print("\nFinal Equation:")
# We display each part of the fraction to show the final calculation clearly.
print(f"({P_all_heads.numerator}/{P_all_heads.denominator}) / ({P_even_tails.numerator}/{P_even_tails.denominator}) = {conditional_probability.numerator}/{conditional_probability.denominator}")
