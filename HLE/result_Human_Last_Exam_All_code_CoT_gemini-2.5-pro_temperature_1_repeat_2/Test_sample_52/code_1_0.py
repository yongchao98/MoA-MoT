from fractions import Fraction

# Step 1: Define the basic probabilities
p_h = Fraction(1, 3)  # Probability of heads
p_t = Fraction(2, 3)  # Probability of tails

# Step 2: Calculate the probability of our event of interest (A: all heads)
# This is also P(A and B) because 0 tails is an even number.
p_all_heads = p_h ** 3

# Step 3: Calculate the probability of our condition (B: even number of tails)
# This can be 0 tails (HHH) or 2 tails (TTH, THT, HTT).
p_0_tails = p_all_heads
# There are 3 ways to get 2 tails, each with probability p_t**2 * p_h
p_2_tails = 3 * (p_t ** 2 * p_h)
p_even_tails = p_0_tails + p_2_tails

# Step 4: Calculate the conditional probability P(A|B) = P(A and B) / P(B)
result = p_all_heads / p_even_tails

# Print the final equation with all the numbers
print("The problem is to find P(All Heads | Even Number of Tails).")
print("This is calculated as P(All Heads and Even Tails) / P(Even Tails).\n")
print("P(All Heads and Even Tails) = P(0 tails) = (1/3)^3 = {}".format(p_all_heads))
print("P(Even Tails) = P(0 tails) + P(2 tails)")
print("P(2 tails) = 3 * (P(Tails))^2 * P(Heads) = 3 * (2/3)^2 * (1/3) = {}".format(p_2_tails))
print("P(Even Tails) = {} + {} = {}".format(p_0_tails, p_2_tails, p_even_tails))
print("\nFinal Calculation:")
print("P(All Heads | Even Tails) = {} / {}".format(p_all_heads, p_even_tails))
print("Result = {}".format(result))

# Final answer in the required format
final_answer_value = float(result)