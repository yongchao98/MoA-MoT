from fractions import Fraction

# Step 1: Define the probabilities for a single coin.
p_h = Fraction(1, 3)
p_t = Fraction(2, 3)

# Step 2: Calculate the probability of the target event (all heads).
# This is also the probability of getting 0 tails.
p_hhh = p_h ** 3

# Step 3: Calculate the probability of getting 2 tails.
# There are 3 combinations for 2 tails: HTT, THT, TTH.
# The probability of any one of these combinations is P(H) * P(T) * P(T).
p_2_tails_one_outcome = p_h * (p_t ** 2)
# The total probability for 2 tails is 3 times this value.
p_2_tails_total = 3 * p_2_tails_one_outcome

# Step 4: Calculate the total probability of the condition (even number of tails).
# This is the sum of the probability of 0 tails (p_hhh) and 2 tails (p_2_tails_total).
p_even_tails = p_hhh + p_2_tails_total

# Step 5: Calculate the final conditional probability.
# P(all heads | even tails) = P(all heads) / P(even tails)
final_probability = p_hhh / p_even_tails

# Step 6: Print the explanation and the final equation.
print("This problem asks for a conditional probability.")
print("Let A be the event 'all three coins are heads' (HHH).")
print("Let B be the event 'the number of tails is even'.")
print("We want to find P(A|B) = P(A and B) / P(B).\n")

print(f"The probability of heads P(H) is {p_h}, and the probability of tails P(T) is {p_t}.\n")

print("First, let's find the probability of event A, P(A).")
print(f"P(A) = P(HHH) = ({p_h}) * ({p_h}) * ({p_h}) = {p_hhh}\n")

print("Next, let's find the probability of event B, P(B).")
print("An even number of tails means 0 tails (HHH) or 2 tails (HTT, THT, TTH).")
print(f"P(0 tails) = P(HHH) = {p_hhh}")
print(f"P(2 tails) = 3 * P(H) * P(T)^2 = 3 * ({p_h}) * ({p_t})^2 = {p_2_tails_total}")
print(f"P(B) = P(0 tails) + P(2 tails) = {p_hhh} + {p_2_tails_total} = {p_even_tails}\n")

print("The event 'A and B' is 'all heads AND an even number of tails'.")
print("This is just the event 'all heads', since 0 is an even number.")
print(f"So, P(A and B) = P(HHH) = {p_hhh}\n")

print("Finally, we can calculate the conditional probability P(A|B).")
print("Final Equation:")
print(f"P(A|B) = P(A and B) / P(B) = {p_hhh} / {p_even_tails} = {final_probability}")
<<<1/13>>>