from fractions import Fraction

# Step 1: Define the probabilities for a single coin toss.
p_h = Fraction(1, 3)  # Probability of heads
p_t = Fraction(2, 3)  # Probability of tails

# Step 2: Calculate the probability of the components of the condition (event B: even number of tails).
# Event B consists of 0 tails or 2 tails.
# Probability of 0 tails (HHH). This is also our event A.
p_0_tails = p_h ** 3

# Probability of 2 tails (e.g., HTT). There are 3 such combinations (HTT, THT, TTH).
p_2_tails = 3 * p_h * (p_t ** 2)

# Total probability of event B (even number of tails).
p_even_tails = p_0_tails + p_2_tails

# Step 3: The event 'A and B' is 'all heads and even tails', which is just 'all heads'.
# So P(A and B) is p_0_tails.
p_a_and_b = p_0_tails

# Step 4: Calculate the final conditional probability P(A|B) = P(A and B) / P(B).
final_probability = p_a_and_b / p_even_tails

# Print the result, showing each number in the final equation as requested.
print("This problem is solved using the conditional probability formula: P(A|B) = P(A and B) / P(B)")
print("A = All three coins are heads.")
print("B = The number of tails is even (0 or 2 tails).\n")
print(f"P(A) = Probability of 0 tails = (1/3)^3 = {p_0_tails}")
print(f"P(B) = P(0 tails) + P(2 tails) = {p_0_tails} + {p_2_tails} = {p_even_tails}\n")
print(f"The required probability P(A|B) is P(A) / P(B).")
print("The final equation is:")
print(f"{p_a_and_b} / {p_even_tails} = {final_probability}")
