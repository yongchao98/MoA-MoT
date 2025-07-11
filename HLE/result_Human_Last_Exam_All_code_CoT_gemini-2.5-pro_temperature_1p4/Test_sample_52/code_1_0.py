from fractions import Fraction

# Step 1: Define the probabilities for a single coin.
# P(H) is the probability of heads.
# P(T) is the probability of tails.
p_h = Fraction(1, 3)
p_t = 1 - p_h

print("This script calculates the conditional probability based on the given problem.\n")
print(f"The probability of a single coin turning heads is P(H) = {p_h}.")
print(f"The probability of a single coin turning tails is P(T) = 1 - P(H) = {p_t}.\n")

# Step 2: Define and calculate the probability of the event of interest (A).
# Event A: All three coins turn heads (HHH).
p_A = p_h ** 3
print("Let A be the event that all three coins turn heads (HHH).")
print(f"The probability of A is P(A) = P(H)^3 = ({p_h})^3 = {p_A}.\n")


# Step 3 & 4: Define and calculate the probability of the conditioning event (E).
# Event E: The number of tails is even (0 or 2 tails).
# Probability of 0 tails (outcome HHH):
p_0_tails = p_h ** 3
# Probability of 2 tails (outcomes TTH, THT, HTT). There are 3 such combinations.
# The probability of one such combination (e.g., TTH) is P(T)*P(T)*P(H).
p_2_tails = 3 * (p_t ** 2) * p_h

# The total probability of event E is the sum.
p_E = p_0_tails + p_2_tails
print("Let E be the event that the number of tails is even.")
print("This occurs if there are 0 tails or 2 tails.")
print(f"The probability of 0 tails (HHH) is {p_0_tails}.")
print(f"The probability of 2 tails (TTH, THT, HTT) is 3 * P(T)^2 * P(H) = 3 * ({p_t})^2 * ({p_h}) = {p_2_tails}.")
print(f"The total probability of E is P(E) = P(0 tails) + P(2 tails) = {p_0_tails} + {p_2_tails} = {p_E}.\n")


# Step 5 & 6: Apply the conditional probability formula.
# We want to find P(A|E) = P(A and E) / P(E).
# Since event A (all heads) is a specific case of event E (0 is an even number of tails),
# the event "A and E" is the same as event A. So, P(A and E) = P(A).
final_probability = p_A / p_E
print("We need to find the probability of A given E, which is P(A|E) = P(A and E) / P(E).")
print(f"Since A is a subset of E, P(A and E) is just P(A).")
print(f"So, P(A|E) = P(A) / P(E).\n")
print("The final calculation is:")
print(f"({p_A.numerator}/{p_A.denominator}) / ({p_E.numerator}/{p_E.denominator}) = {final_probability.numerator}/{final_probability.denominator}")
print("\nIn full, the final equation with each number is:")
print(f"P(A|E) = P(HHH) / (P(HHH) + P(2 Tails)) = ({p_A.numerator}/{p_A.denominator}) / ({p_0_tails.numerator}/{p_0_tails.denominator} + {p_2_tails.numerator}/{p_2_tails.denominator}) = {final_probability.numerator}/{final_probability.denominator}")

<<<1/13>>>