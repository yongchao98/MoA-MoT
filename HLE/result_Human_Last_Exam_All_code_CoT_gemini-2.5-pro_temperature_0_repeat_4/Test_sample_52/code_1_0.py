from fractions import Fraction

# Step 1: Define the base probabilities for a single coin.
p_h = Fraction(1, 3)
p_t = 1 - p_h

# Step 2: Calculate the probability of the target event (all heads).
# This is also the event of having 0 tails.
p_all_heads = p_h ** 3

# Step 3: Calculate the probability of the condition (even number of tails).
# This can be 0 tails (HHH) or 2 tails (HTT, THT, TTH).

# The probability of 2 tails and 1 head can happen in 3 ways (HTT, THT, TTH).
# The probability of one such specific outcome (e.g., HTT) is:
p_one_way_for_2_tails = p_h * p_t * p_t
# The total probability for 2 tails is 3 times this value.
p_total_for_2_tails = 3 * p_one_way_for_2_tails

# The total probability of the condition (even tails) is the sum of
# the probability of 0 tails and the probability of 2 tails.
p_condition_even_tails = p_all_heads + p_total_for_2_tails

# Step 4: Apply the conditional probability formula.
# P(All Heads | Even Tails) = P(All Heads) / P(Even Tails)
final_probability = p_all_heads / p_condition_even_tails

# Step 5: Print the step-by-step calculation.
print("This script calculates the probability of getting all heads given that the number of tails is even.")
print("-" * 80)
print(f"The probability of a single head P(H) is {p_h}.")
print(f"The probability of a single tail P(T) is {p_t}.")
print("-" * 80)
print("First, let's calculate the probability of the event 'all heads' (HHH).")
print(f"P(HHH) = P(H) * P(H) * P(H) = {p_h} * {p_h} * {p_h} = {p_all_heads}")
print("-" * 80)
print("Next, let's calculate the probability of the condition 'the number of tails is even'.")
print("This occurs with 0 tails (HHH) or 2 tails (HTT, THT, TTH).")
print(f"P(0 tails) is the same as P(HHH), which is {p_all_heads}.")
print(f"P(2 tails) = 3 * P(H) * P(T) * P(T) = 3 * {p_h} * {p_t} * {p_t} = {p_total_for_2_tails}")
print(f"P(Even Tails) = P(0 tails) + P(2 tails) = {p_all_heads} + {p_total_for_2_tails} = {p_condition_even_tails}")
print("-" * 80)
print("Finally, we calculate the conditional probability P(All Heads | Even Tails).")
print("The formula is P(HHH) / P(Even Tails).")
print(f"Final Equation: {p_all_heads} / {p_condition_even_tails}")
print(f"The final probability is: {final_probability}")
<<<1/13>>>