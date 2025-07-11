import numpy as np
from fractions import Fraction

def solve_bran_castle_escape():
    """
    Calculates the probability of reaching the Treasure Room from the Main Entrance
    using absorbing Markov chain analysis.
    """
    
    # The states are ordered as follows:
    # 0: Main Entrance (ME)
    # 1: Great Hall (GH)
    # 2: Castle Kitchen (CK)
    # 3: Secret Passage (SP)
    # 4: Library (L)
    # 5: Knights' Hall (KH)
    # 6: Torture Chamber (TC)
    # 7: Queen's Room (QR)
    #
    # The absorbing states are:
    # 0: Treasure Room (TR)
    # 1: Vampire's Lair (VL)

    # Q is the transition matrix for transient states (transient -> transient)
    Q = np.array([
        # ME   GH        CK        SP        L         KH        TC        QR
        [0, Fraction(1,2), Fraction(1,2), 0,        0,        0,        0,        0],         # ME
        [Fraction(1,4), 0,        0,        Fraction(1,4), Fraction(1,4), Fraction(1,4), 0,        0],         # GH
        [Fraction(1,2), 0,        0,        0,        0,        Fraction(1,2), 0,        0],         # CK
        [0, Fraction(1,4), 0,        0,        0,        0,        Fraction(1,4), 0],         # SP
        [0, Fraction(1,3), 0,        0,        0,        Fraction(1,3), 0,        Fraction(1,3)], # L
        [0, Fraction(1,4), Fraction(1,4), 0,        Fraction(1,4), 0,        0,        Fraction(1,4)], # KH
        [0, 0,        0,        Fraction(1,2), 0,        0,        0,        0],         # TC
        [0, 0,        0,        0,        Fraction(1,3), Fraction(1,3), 0,        0]          # QR
    ], dtype=object)

    # R is the transition matrix from transient to absorbing states
    R = np.array([
        # TR        VL
        [0,        0],           # ME
        [0,        0],           # GH
        [0,        0],           # CK
        [Fraction(1,4), Fraction(1,4)], # SP
        [0,        0],           # L
        [0,        0],           # KH
        [0,        Fraction(1,2)], # TC
        [Fraction(1,3), 0]            # QR
    ], dtype=object)

    # The fundamental matrix N = (I - Q)^-1
    I = np.identity(len(Q), dtype=object)
    I_minus_Q = I - Q
    N = np.linalg.inv(I_minus_Q)

    # The absorption probability matrix B = N * R
    B = N @ R

    # The probability of reaching the Treasure Room (absorbing state 0)
    # starting from the Main Entrance (transient state 0) is B[0, 0].
    probability = B[0, 0]
    
    print("This problem can be modeled by a system of linear equations, where P_room is the probability of success from that room.")
    print("The system of equations is:")
    print("P_ME = (1/2) * P_GH + (1/2) * P_CK")
    print("P_GH = (1/4) * P_ME + (1/4) * P_SP + (1/4) * P_L + (1/4) * P_KH")
    print("P_CK = (1/2) * P_ME + (1/2) * P_KH")
    print("P_SP = (1/4) * P_GH + (1/4) * P_TC + 1/4")
    print("P_L = (1/3) * P_GH + (1/3) * P_KH + (1/3) * P_QR")
    print("P_KH = (1/4) * P_GH + (1/4) * P_CK + (1/4) * P_L + (1/4) * P_QR")
    print("P_TC = (1/2) * P_SP")
    print("P_QR = (1/3) * P_L + (1/3) * P_KH + 1/3")
    print("\nSolving this system for P_ME gives the final probability.")
    
    # The final equation is the value of P_ME
    print("\nThe probability of reaching the Treasure Room from the Main Entrance is:")
    print(f"P_ME = {probability.numerator} / {probability.denominator}")


solve_bran_castle_escape()