from fractions import Fraction

# Step 1 & 2: Define probabilities
p_h = Fraction(1, 3)
p_t = Fraction(2, 3)

# Step 3: Calculate the probability of the events

# Event A: All coins are heads (HHH). This is also the case for 0 tails.
p_all_heads = p_h**3

# Event B: The number of tails is even (0 or 2 tails).
# Probability of 0 tails is p_all_heads.
p_0_tails = p_all_heads

# Probability of 2 tails. There are 3 combinations (HTT, THT, TTH).
# Probability of one combination, e.g., HTT, is P(H)*P(T)*P(T).
p_one_combo_2_tails = p_h * (p_t**2)
p_2_tails = 3 * p_one_combo_2_tails

# Total probability of an even number of tails
p_even_tails = p_0_tails + p_2_tails

# Step 4: Calculate the final conditional probability.
# P(all heads | even tails) = P(all heads) / P(even tails)
final_probability = p_all_heads / p_even_tails

# Step 5: Print the explanation and the result.
print("This problem asks for the conditional probability P(all heads | number of tails is even).\n")
print(f"First, we find the probability of all heads (0 tails):")
print(f"P(all heads) = (1/3) * (1/3) * (1/3) = {p_all_heads.numerator}/{p_all_heads.denominator}\n")

print(f"Next, we find the probability of the number of tails being even (0 or 2 tails):")
print(f"P(0 tails) = {p_0_tails.numerator}/{p_0_tails.denominator}")
print(f"P(2 tails) = 3 * (1/3 * 2/3 * 2/3) = {p_2_tails.numerator}/{p_2_tails.denominator}")
print(f"P(even tails) = P(0 tails) + P(2 tails) = {p_0_tails.numerator}/{p_0_tails.denominator} + {p_2_tails.numerator}/{p_2_tails.denominator} = {p_even_tails.numerator}/{p_even_tails.denominator}\n")

print("Finally, we compute the conditional probability:")
print("P(all heads | even tails) = P(all heads) / P(even tails)")
print(f"The equation is: ({p_all_heads.numerator}/{p_all_heads.denominator}) / ({p_even_tails.numerator}/{p_even_tails.denominator}) = {final_probability.numerator}/{final_probability.denominator}")
print(f"\nThe result is {final_probability.numerator}/{final_probability.denominator}, which is approximately {float(final_probability):.4f}.")
