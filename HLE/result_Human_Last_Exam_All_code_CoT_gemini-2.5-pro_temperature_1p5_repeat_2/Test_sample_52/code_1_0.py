import math

# Define probabilities for a single coin
p_h = 1/3
p_t = 1 - p_h  # Probability of tails is 1 - P(Heads)

# Number of coins
n = 3

# --- Calculation for Event A (All Heads) ---
# P(A) = P(HHH) = P(H)^3
# For printing, we calculate the numerator and denominator separately.
prob_A_numerator = 1**n
prob_A_denominator = 3**n
prob_A = prob_A_numerator / prob_A_denominator

# --- Calculation for Event B (Even number of tails: 0 or 2) ---
# P(B) = P(0 tails) + P(2 tails)

# P(0 tails) is the same as P(A)
prob_0_tails_num = prob_A_numerator
prob_0_tails_den = prob_A_denominator

# P(2 tails) = C(3, 2) * P(H)^1 * P(T)^2
# C(n, k) is the number of combinations 'n choose k'
num_combinations_2_tails = math.comb(n, 2)
# For P(TTH) = P(T)^2 * P(H)^1, the numerator is 2^2 * 1^1 = 4
prob_one_combination_2_tails_num = (2**2) * (1**1)
prob_2_tails_num = num_combinations_2_tails * prob_one_combination_2_tails_num
prob_2_tails_den = 3**n

# P(B) numerator = P(0 tails) numerator + P(2 tails) numerator
prob_B_numerator = prob_0_tails_num + prob_2_tails_num
prob_B_denominator = prob_A_denominator
prob_B = prob_B_numerator / prob_B_denominator

# --- Final Conditional Probability P(A|B) ---
# P(A|B) = P(A) / P(B)
# The denominators (27) cancel out.
final_prob_numerator = prob_A_numerator
final_prob_denominator = prob_B_numerator
final_probability = prob_A / prob_B

# --- Print the explanation and final equation ---
print("Let A be the event of all heads (HHH) and B be the event of an even number of tails.")
print("We need to find the conditional probability P(A|B) = P(A and B) / P(B).")
print("Since getting all heads means 0 tails (an even number), event A implies event B.")
print("Therefore, P(A and B) = P(A).\n")

print("1. Calculate the probability of event A (all heads):")
print(f"P(A) = P(H)^3 = (1/3)^3 = {prob_A_numerator}/{prob_A_denominator}\n")

print("2. Calculate the probability of event B (even number of tails):")
print("P(B) = P(0 tails) + P(2 tails)")
print(f"P(0 tails) = P(A) = {prob_0_tails_num}/{prob_0_tails_den}")
print(f"P(2 tails) = C(3, 2) * P(H)^1 * P(T)^2 = 3 * (1/3) * (2/3)^2 = {prob_2_tails_num}/{prob_2_tails_den}")
print(f"P(B) = {prob_0_tails_num}/{prob_0_tails_den} + {prob_2_tails_num}/{prob_2_tails_den} = {prob_B_numerator}/{prob_B_denominator}\n")

print("3. Calculate the final probability P(A|B):")
print("P(A|B) = P(A) / P(B)")
print(f"       = ({prob_A_numerator}/{prob_A_denominator}) / ({prob_B_numerator}/{prob_B_denominator})")
print(f"       = {final_prob_numerator}/{final_prob_denominator}")

<<<1/13>>>