from fractions import Fraction

# Define the probability for a single coin flip
p_head = Fraction(1, 3)
p_tail = Fraction(2, 3)

# 1. Calculate the probability of the event 'all heads' (HHH), which corresponds to 0 tails.
# P(HHH) = P(H) * P(H) * P(H)
p_all_heads = p_head**3

# 2. Calculate the probability of getting 2 tails.
# There are 3 combinations for 2 tails and 1 head (HTT, THT, TTH).
# The probability for any one of these specific combinations is P(H) * P(T) * P(T).
p_one_combo_2_tails = p_head * (p_tail**2)
# The total probability of getting 2 tails is 3 times this value.
p_total_2_tails = 3 * p_one_combo_2_tails

# 3. Calculate the total probability of the condition "number of tails is even".
# This is P(0 tails) + P(2 tails) = P(HHH) + P(2 tails).
p_even_tails = p_all_heads + p_total_2_tails

# 4. Calculate the final conditional probability.
# P(All Heads | Even Tails) = P(All Heads) / P(Even Tails)
final_probability = p_all_heads / p_even_tails

print("The final probability is calculated using the formula:")
print("P(All Heads | Even Tails) = P(All Heads) / (P(0 Tails) + P(2 Tails))")
print("\nSubstituting the calculated values into the equation:")
print(f"P(All Heads | Even Tails) = ({p_all_heads}) / ({p_all_heads} + {p_total_2_tails})")
print(f"                             = ({p_all_heads}) / ({p_even_tails})")
print(f"                             = {final_probability}")

# <<<1/13>>>