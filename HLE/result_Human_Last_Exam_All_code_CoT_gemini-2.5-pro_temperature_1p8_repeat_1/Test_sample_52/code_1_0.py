import math

def combinations(n, k):
    """Calculates the number of combinations 'n choose k'."""
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

# Probabilities for a single coin
p_h = 1/3
p_t = 2/3

# Number of coins
n = 3

print("Let's calculate the required probabilities:\n")

# --- Probability of Event A and B: All Heads (0 Tails) ---
# Event A is 'all heads' (HHH). This means there are 0 tails, which is an even number.
# So, the event 'A and B' (all heads AND an even number of tails) is just the same as event A.
# P(A and B) = P(HHH) = P(H) * P(H) * P(H)
p_all_heads_num = int(p_h.denominator) ** n
p_all_heads = 1 / p_all_heads_num
print(f"1. Probability of all heads (0 tails):")
print(f"   P(HHH) = (1/3)^3 = 1/{p_all_heads_num}")
print("-" * 30)

# --- Probability of Event B: Even Number of Tails (0 or 2) ---
# The number of tails can be even in two ways: 0 tails or 2 tails.

# We already have P(0 tails) = 1/27.
p_0_tails = p_all_heads

# Now, let's calculate the probability of getting exactly 2 tails.
# The number of ways to get 2 tails from 3 coins is C(3, 2).
# Each specific combination (e.g., HTT) has a probability of P(H)^1 * P(T)^2.
num_ways_2_tails = combinations(n, n - 1)
p_one_of_2_tails = p_h * (p_t ** 2)
p_2_tails = num_ways_2_tails * p_one_of_2_tails

# Get numerators and denominators for printing
p_2_tails_num = num_ways_2_tails * int(p_t.numerator ** 2)
p_2_tails_den = p_all_heads_num

print(f"2. Probability of an even number of tails:")
print(f"   This means either 0 tails or 2 tails.\n")
print(f"   a) Probability of 2 tails:")
print(f"      - Number of combinations for 2 tails (e.g., HTT, THT, TTH) = C(3, 2) = {num_ways_2_tails}")
print(f"      - Probability of one such combination = (1/3) * (2/3)^2 = 4/27")
print(f"      - Total P(2 tails) = {num_ways_2_tails} * (4/27) = {p_2_tails_num}/{p_2_tails_den}\n")

# P(B) = P(0 tails) + P(2 tails)
p_even_tails = p_0_tails + p_2_tails
p_even_tails_num = 1 + p_2_tails_num
p_even_tails_den = p_all_heads_num

print(f"   b) Total probability of even tails:")
print(f"      P(even tails) = P(0 tails) + P(2 tails)")
print(f"      P(even tails) = 1/{p_all_heads_num} + {p_2_tails_num}/{p_2_tails_den} = {p_even_tails_num}/{p_even_tails_den}")
print("-" * 30)


# --- Final Conditional Probability ---
# P(A|B) = P(A and B) / P(B)
final_prob = p_all_heads / p_even_tails
final_prob_num = 1
final_prob_den = p_even_tails_num

print("3. Final conditional probability:")
print("   P(all heads | even tails) = P(all heads) / P(even tails)\n")
print(f"   The final equation is:")
print(f"   (1/{p_all_heads_num}) / ({p_even_tails_num}/{p_even_tails_den}) = {final_prob_num}/{final_prob_den}\n")
print(f"The final probability is {final_prob_num}/{final_prob_den}, which is approximately {final_prob:.4f}")

<<<1/13>>>