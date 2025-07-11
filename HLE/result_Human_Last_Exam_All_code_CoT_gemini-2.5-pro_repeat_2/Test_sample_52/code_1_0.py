import fractions

# Step 1: Define base probabilities
p_h = fractions.Fraction(1, 3)
p_t = 1 - p_h

# The sample space is restricted to outcomes with an even number of tails (0 or 2).

# Step 2: Calculate the probability of the event of interest: 0 tails (all heads).
# This corresponds to the outcome HHH.
p_0_tails = p_h ** 3

# Step 3: Calculate the probability of the other part of the condition: 2 tails.
# This corresponds to the outcomes HTT, THT, TTH. There are 3 such outcomes.
# The probability of any single one of these (e.g., HTT) is p_h * p_t * p_t.
p_2_tails = 3 * p_h * (p_t ** 2)

# Step 4: Calculate the total probability of the condition (even number of tails).
# This is the sum of the probabilities from steps 2 and 3.
p_even_tails = p_0_tails + p_2_tails

# Step 5: Calculate the conditional probability.
# P(all heads | even tails) = P(all heads AND even tails) / P(even tails)
# Since "all heads" (0 tails) is an event with an even number of tails,
# P(all heads AND even tails) is simply P(all heads).
final_probability = p_0_tails / p_even_tails

# Print the final equation with all the numbers calculated.
print("The probability of all heads given an even number of tails is calculated as:")
print("P(all heads) / P(even tails)")
print("\nWhere:")
print(f"P(all heads) = P(0 tails) = (1/3)^3 = {p_0_tails}")
print(f"P(2 tails) = 3 * (1/3) * (2/3)^2 = {p_2_tails}")
print(f"P(even tails) = P(0 tails) + P(2 tails) = {p_0_tails} + {p_2_tails} = {p_even_tails}")
print("\nSo the final equation is:")
print(f"{p_0_tails} / {p_even_tails} = {final_probability}")
print(f"\nThe final probability is {final_probability.numerator}/{final_probability.denominator}")
