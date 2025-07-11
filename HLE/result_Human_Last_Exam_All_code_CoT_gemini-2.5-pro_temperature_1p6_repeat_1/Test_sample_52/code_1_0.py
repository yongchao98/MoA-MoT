import fractions

# Step 1: Define the initial probabilities
p_h = fractions.Fraction(1, 3)
p_t = 1 - p_h

# Step 2: Calculate the probability of the outcomes that satisfy the condition (even tails)

# Case 1: 0 tails (HHH)
p_0_tails = p_h * p_h * p_h

# Case 2: 2 tails (HTT, THT, TTH)
# The number of ways to get 2 tails in 3 tosses is 3.
num_ways_2_tails = 3
p_2_tails_one_way = p_h * p_t * p_t
p_2_tails_total = num_ways_2_tails * p_2_tails_one_way

# Step 3: Calculate the total probability of the condition (even number of tails)
# This is our new sample space probability.
p_even_tails = p_0_tails + p_2_tails_total

# Step 4: Calculate the conditional probability
# The event of "all heads" (HHH) is the same as "0 tails".
# P(all heads | even tails) = P(all heads AND even tails) / P(even tails)
# Since "all heads" has 0 tails (which is even), P(all heads AND even tails) is just P(all heads).
# So, P(all heads | even tails) = P(0 tails) / P(even tails)
final_prob = p_0_tails / p_even_tails

# Print the explanation and the result
print("Let P(H) be the probability of heads and P(T) be the probability of tails.")
print(f"P(H) = {p_h}")
print(f"P(T) = {p_t}")
print("\nThe condition is that the number of tails is even (0 or 2 tails).")
print("\n1. Probability of 0 tails (HHH):")
print(f"P(HHH) = ({p_h}) * ({p_h}) * ({p_h}) = {p_0_tails}")

print("\n2. Probability of 2 tails (HTT, THT, TTH):")
print(f"There are 3 ways to get 2 tails. The probability of one way (e.g., HTT) is ({p_h})*({p_t})*({p_t}) = {p_2_tails_one_way}.")
print(f"Total probability for 2 tails = 3 * {p_2_tails_one_way} = {p_2_tails_total}")

print("\n3. Total probability of the number of tails being even:")
print(f"P(even tails) = P(0 tails) + P(2 tails) = {p_0_tails} + {p_2_tails_total} = {p_even_tails}")

print("\n4. The desired conditional probability is P(all heads | even tails):")
print("P(all heads | even tails) = P(all heads) / P(even tails)")

# In the final equation, we show the numbers that lead to the final fraction
print(f"\nFinal Equation: {p_0_tails.numerator}/{p_0_tails.denominator} / ({p_even_tails.numerator}/{p_even_tails.denominator}) = {final_prob.numerator} / {final_prob.denominator}")
print(f"The final probability is {final_prob}.")
