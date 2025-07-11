import math

# Define the probabilities for a single coin toss
P_H = 1/3
P_T = 1 - P_H

# --- 1. Calculate the numerator: P(A and B) ---
# Event A = All heads (HHH)
# Event B = Even number of tails (0 or 2)
# Event (A and B) = HHH, because it has 0 tails (which is even).
# We calculate its probability as a fraction to maintain precision.
prob_A_and_B_num = 1**3
prob_A_and_B_den = 3**3

print(f"The probability of all heads (HHH), which has 0 tails, is P(A and B).")
print(f"P(A and B) = P(H)^3 = (1/3)^3 = {prob_A_and_B_num}/{prob_A_and_B_den}\n")

# --- 2. Calculate the denominator: P(B) ---
# P(B) = P(0 tails) + P(2 tails)
# P(0 tails) is the same as P(A and B)
prob_0_tails_num = prob_A_and_B_num
prob_0_tails_den = prob_A_and_B_den

# P(2 tails) = C(3, 2) * P(H)^1 * P(T)^2
# C(3, 2) is the number of ways to get 2 tails from 3 coins.
num_ways_2_tails = math.comb(3, 2)
prob_2_tails_num = num_ways_2_tails * (1**1) * (2**2)
prob_2_tails_den = 3**3

# P(B) is the sum of these probabilities
prob_B_num = prob_0_tails_num + prob_2_tails_num
prob_B_den = prob_2_tails_den # Denominators are the same

print(f"The probability of an even number of tails, P(B), is P(0 tails) + P(2 tails).")
print(f"P(0 tails) = {prob_0_tails_num}/{prob_0_tails_den}")
print(f"P(2 tails) = C(3,2) * P(H)^1 * P(T)^2 = {num_ways_2_tails} * (1/3) * (4/9) = {prob_2_tails_num}/{prob_2_tails_den}")
print(f"P(B) = {prob_0_tails_num}/{prob_0_tails_den} + {prob_2_tails_num}/{prob_2_tails_den} = {prob_B_num}/{prob_B_den}\n")


# --- 3. Calculate the final conditional probability P(A|B) ---
# P(A|B) = P(A and B) / P(B)
final_num = prob_A_and_B_num
final_den = prob_B_num

# Simplify the fraction using the greatest common divisor if needed, but in this case, the numerator is 1.
common_divisor = math.gcd(final_num, final_den)
simplified_num = final_num // common_divisor
simplified_den = final_den // common_divisor


print("The final probability of all heads given an even number of tails is:")
print(f"P(A|B) = P(A and B) / P(B) = ({prob_A_and_B_num}/{prob_A_and_B_den}) / ({prob_B_num}/{prob_B_den})")
print(f"P(A|B) = {prob_A_and_B_num}/{prob_B_num} = {simplified_num}/{simplified_den}")
<<<1/13>>>