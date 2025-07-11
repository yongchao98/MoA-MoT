from fractions import Fraction

# Step 1: Define the base probabilities for a single coin.
# Using the Fraction class to avoid floating-point inaccuracies and show the math clearly.
p_h = Fraction(1, 3)
p_t = 1 - p_h

# Step 2: Define the event of interest and the conditional event.
# The event we want the probability of is "all heads" (HHH).
# This outcome also satisfies the condition (0 tails is even), so it's our "A and B" event.
p_event_of_interest = p_h ** 3

# The condition is that the number of tails is even (0 or 2 tails).
# Calculate the probability of 0 tails (which is the same as our event of interest).
p_0_tails = p_event_of_interest

# Calculate the probability of 2 tails. There are 3 ways this can happen: HTT, THT, TTH.
# The probability of any one of these specific outcomes is p_h * p_t * p_t.
# Number of combinations for 2 tails out of 3 tosses is 3.
num_combinations_2_tails = 3
p_2_tails = num_combinations_2_tails * p_h * (p_t ** 2)

# Step 3: Calculate the total probability of the condition (our new sample space).
p_condition = p_0_tails + p_2_tails

# Step 4: Calculate the final conditional probability.
# P(A|B) = P(A and B) / P(B)
# Here, A is "all heads" and B is "even number of tails".
# P(A and B) is p_event_of_interest. P(B) is p_condition.
final_probability = p_event_of_interest / p_condition

# Step 5: Print the explanation and the final equation with all numbers.
print("This is a conditional probability problem.")
print("Let A be the event 'all heads' (HHH).")
print("Let B be the event 'the number of tails is even' (0 or 2 tails).\n")

print("We want to find P(A|B) = P(A and B) / P(B).\n")

print("The probability of event A (HHH), which has 0 tails (an even number), is:")
print(f"P(A and B) = P(HHH) = (1/3)^3 = {p_event_of_interest}\n")

print("The probability of event B (even tails) is P(0 tails) + P(2 tails).")
print(f"P(2 tails) = 3 * P(one H, two T) = 3 * (1/3) * (2/3)^2 = {p_2_tails}")
print(f"P(B) = P(0 tails) + P(2 tails) = {p_0_tails} + {p_2_tails} = {p_condition}\n")

print("The final calculation is:")
print(f"P(A|B) = P(A and B) / P(B) = {p_event_of_interest} / {p_condition} = {final_probability}")
<<<1/13>>>