from fractions import Fraction

# Step 1: Define the probabilities for a single coin toss.
# P(H) is the probability of a head.
# P(T) is the probability of a tail.
p_h = Fraction(1, 3)
p_t = 1 - p_h

# Step 2: Calculate the probability of the condition B: "the number of tails is even".
# This means 0 tails or 2 tails.

# Probability of 0 tails (HHH)
p_0_tails = p_h ** 3

# Probability of 2 tails (HTT, THT, TTH).
# There are 3 such outcomes, each with probability P(H) * P(T)^2.
p_2_tails = 3 * p_h * (p_t ** 2)

# Total probability of the condition B (even number of tails).
p_even_tails = p_0_tails + p_2_tails

# Step 3: Calculate the probability of the event A: "all three coins are heads".
# This is the outcome HHH, which is the same as P(0 tails).
p_all_heads = p_0_tails

# Step 4: The problem asks for P(A | B), the probability of A given B.
# The formula is P(A | B) = P(A and B) / P(B).
# The event "A and B" is "all heads AND an even number of tails".
# Since "all heads" (HHH) has 0 tails, which is an even number, the event A is a subset of B.
# Therefore, P(A and B) is the same as P(A).
p_intersection = p_all_heads

# Step 5: Calculate the final conditional probability.
conditional_probability = p_intersection / p_even_tails

# Step 6: Print the full calculation.
print("This problem is about conditional probability.")
print("We want to find P(All Heads | Even Number of Tails).\n")
print("Let P(H) = 1/3 and P(T) = 2/3.\n")

print("First, we find the probability of the condition (even number of tails):")
print(f"P(0 Tails) = P(HHH) = (1/3)^3 = {p_0_tails}")
print(f"P(2 Tails) = 3 * P(HTT) = 3 * (1/3) * (2/3)^2 = {p_2_tails}")
print(f"P(Even Tails) = P(0 Tails) + P(2 Tails) = {p_0_tails} + {p_2_tails} = {p_even_tails}\n")

print("The event 'All Heads' (HHH) has 0 tails, which is an even number.")
print("So, P(All Heads and Even Tails) = P(All Heads) = P(HHH).\n")

print("Finally, we calculate the conditional probability:")
print("P(All Heads | Even Tails) = P(All Heads and Even Tails) / P(Even Tails)")

# We use .numerator and .denominator to clearly show the fractions in the final equation.
print(f"The final equation is: ({p_intersection.numerator}/{p_intersection.denominator}) / ({p_even_tails.numerator}/{p_even_tails.denominator}) = {conditional_probability.numerator}/{conditional_probability.denominator}")

print(f"\nThe probability is {conditional_probability}.")

<<<1/13>>>