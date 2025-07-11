from fractions import Fraction

# Step 1: Define the base probabilities for a single coin toss.
P_H = Fraction(1, 3)
P_T = 1 - P_H  # This is 2/3

# The problem asks for the conditional probability:
# P(all heads | number of tails is even)

# This can be written as P(A | B) = P(A and B) / P(B)
# Event A = All coins are heads (HHH).
# Event B = The number of tails is even (i.e., 0 tails or 2 tails).

# Step 2: Calculate P(A and B).
# The event "all heads" (A) means there are 0 tails. Since 0 is an even number,
# this event satisfies condition B. Therefore, the event "A and B" is just A.
# P(A and B) = P(A) = P(HHH)
p_A_and_B = P_H ** 3

# Step 3: Calculate P(B), the probability of the condition.
# P(B) = P(0 tails) + P(2 tails)

# P(0 tails) is the same as P(all heads).
p_0_tails = P_H ** 3

# P(2 tails) can occur in 3 ways: HTT, THT, TTH.
# Each combination has one Head and two Tails.
# The number of combinations is C(3, 2) = 3.
p_2_tails = 3 * (P_H ** 1) * (P_T ** 2)

# Now, add them to get P(B).
p_B = p_0_tails + p_2_tails

# Step 4: Calculate the final conditional probability P(A|B).
conditional_prob = p_A_and_B / p_B

# Step 5: Print the detailed explanation and the final equation.
print("This is a conditional probability problem: P(all heads | even number of tails).")
print("Let P(H) = 1/3 and P(T) = 2/3.\n")

print("The probability of the event we want ('all heads') is:")
print(f"P(all heads) = P(HHH) = (1/3)^3 = {p_A_and_B}\n")

print("The probability of the condition ('even number of tails') is P(0 tails) + P(2 tails).")
print(f"P(0 tails) = P(HHH) = {p_0_tails}")
print(f"P(2 tails) = 3 * P(H) * P(T)^2 = 3 * (1/3) * (2/3)^2 = {p_2_tails}")
print(f"P(even number of tails) = {p_0_tails} + {p_2_tails} = {p_B}\n")

print("The final probability is P(all heads) / P(even number of tails).")
print("The final equation is:")
# The requested output format showing each number in the equation.
print(f"{p_A_and_B} / {p_B} = {conditional_prob}")
