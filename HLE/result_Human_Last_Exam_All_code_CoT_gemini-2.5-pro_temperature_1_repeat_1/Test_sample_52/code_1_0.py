import fractions

# --- Step 1: Define initial probabilities ---
# Use the fractions module for precision
p_h = fractions.Fraction(1, 3)
p_t = 1 - p_h

# --- Step 2: Calculate the probability of the condition (Event B: even number of tails) ---
# Case 1: 0 tails (HHH)
# The probability is P(H) * P(H) * P(H)
prob_0_tails = p_h ** 3

# Case 2: 2 tails (TTH, THT, HTT)
# The probability of one combination (e.g., TTH) is P(T)*P(T)*P(H)
# There are 3 such combinations.
prob_2_tails = 3 * (p_t ** 2 * p_h)

# The total probability of Event B is the sum of these cases
prob_b = prob_0_tails + prob_2_tails

# --- Step 3: Calculate the probability of the intersection (Event A and B) ---
# Event A is all heads (HHH), which means 0 tails.
# Since 0 is an even number, Event A is a subset of Event B.
# Therefore, P(A and B) is the same as P(A).
prob_a_and_b = prob_0_tails

# --- Step 4: Calculate the final conditional probability P(A|B) ---
# P(A|B) = P(A and B) / P(B)
final_prob = prob_a_and_b / prob_b

# --- Step 5: Print the explanation and the final equation ---
print("This is a conditional probability problem: P(All Heads | Even Tails).")
print("The formula is P(A|B) = P(A and B) / P(B)\n")

print(f"The probability of the condition (even tails) is P(B).")
print(f"P(B) = P(0 tails) + P(2 tails) = {prob_0_tails} + {prob_2_tails} = {prob_b}\n")

print(f"The probability of the intersection (all heads AND even tails) is P(A and B).")
print(f"This is the same as the probability of all heads, which is {prob_a_and_b}\n")

print("The final probability is calculated as follows:")
print(f"P(All Heads | Even Tails) = {prob_a_and_b} / {prob_b} = {final_prob}")

<<<1/13>>>