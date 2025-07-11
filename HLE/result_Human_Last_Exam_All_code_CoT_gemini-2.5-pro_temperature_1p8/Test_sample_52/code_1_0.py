import fractions

# Step 1: Define the base probabilities for a single coin toss.
# P_H is the probability of getting a Head.
# P_T is the probability of getting a Tail.
P_H = fractions.Fraction(1, 3)
P_T = 1 - P_H

print("The probability of one coin being Heads, P(H) = 1/3")
print(f"The probability of one coin being Tails, P(T) = 1 - 1/3 = {P_T}\n")

# The problem states that the number of tails is always even. For 3 coins,
# this means we can have either 0 tails or 2 tails.

# Step 2: Calculate the probability of the outcome we're interested in: all heads (0 tails).
# This is our numerator.
prob_all_heads = P_H * P_H * P_H
print("The event 'all heads' (HHH) has 0 tails, which is an even number.")
print(f"The probability of 'all heads' is P(HHH) = P(H)*P(H)*P(H) = {P_H} * {P_H} * {P_H} = {prob_all_heads}\n")

# Step 3: Calculate the probabilities of other outcomes that also have an even number of tails (2 tails).
# There are 3 ways to get 2 tails: HTT, THT, TTH.
prob_one_outcome_with_2_tails = P_H * P_T * P_T
# The total probability of getting 2 tails is 3 times this, as HTT, THT, and TTH are mutually exclusive.
prob_total_2_tails = 3 * prob_one_outcome_with_2_tails
print("The other way to get an even number of tails is to get exactly 2 tails (e.g., HTT, THT, TTH).")
print(f"The probability of any single one of these outcomes (like HTT) is P(H)*P(T)*P(T) = {P_H} * {P_T} * {P_T} = {prob_one_outcome_with_2_tails}")
print(f"The total probability of getting 2 tails is 3 * {prob_one_outcome_with_2_tails} = {prob_total_2_tails}\n")


# Step 4: Calculate the total probability of the condition "the number of tails is even".
# This is the sum of the probabilities of all outcomes that satisfy the condition.
# This will be the denominator in our final calculation.
prob_condition_even_tails = prob_all_heads + prob_total_2_tails
print("The total probability of our known condition ('number of tails is even') is the sum of these possibilities:")
print(f"P(even tails) = P(0 tails) + P(2 tails) = {prob_all_heads} + {prob_total_2_tails} = {prob_condition_even_tails}\n")


# Step 5: Calculate the final conditional probability.
# P(all heads | number of tails is even) = P(all heads) / P(number of tails is even)
final_probability = prob_all_heads / prob_condition_even_tails
print("The question asks for P(all heads | even tails), which is calculated as P(all heads) / P(even tails).")
print("The final equation is:")
print(f"{prob_all_heads} / {prob_condition_even_tails} = {final_probability}\n")

print(f"Therefore, the final probability of them all turning heads given the number of tails is even is {final_probability}.")
<<<1/13>>>