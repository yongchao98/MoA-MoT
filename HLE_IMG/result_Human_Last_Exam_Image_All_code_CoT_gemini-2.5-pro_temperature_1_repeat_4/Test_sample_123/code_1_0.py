import numpy as np
from fractions import Fraction

def solve_bran_castle_challenge():
    """
    This function sets up and solves the system of linear equations for the
    Bran Castle escape challenge.
    """

    # The problem is to find the probability of reaching the Treasure Room.
    # Let P_i be this probability starting from room i.
    # The non-terminal rooms are ordered as follows for our system of equations:
    # 0: Main Entrance (ME)
    # 1: Great Hall (GH)
    # 2: Castle Kitchen (CK)
    # 3: Secret Passage (SP)
    # 4: Library (L)
    # 5: Knights' Hall (KH)
    # 6: Torture Chamber (TC)
    # 7: Queen's Room (QR)

    # For the absorbing states, we have:
    # P_TreasureRoom = 1
    # P_VampiresLair = 0

    # We set up a system of linear equations Ax = b, where x is the vector of probabilities.
    # Each row in the system corresponds to an equation for one room, derived from:
    # P_i = sum(Prob(move to j from i) * P_j) over all adjacent rooms j.

    # For example, for the Main Entrance: P_ME = 1/2 * P_GH + 1/2 * P_CK
    # Rearranged for the matrix form: 1*P_ME - 1/2*P_GH - 1/2*P_CK = 0

    # Matrix A contains the coefficients of the P_i variables.
    # Vector b contains the constant terms, which come from paths leading
    # directly to the Treasure Room or Vampire's Lair.
    
    #       ME,           GH,            CK,            SP,            L,             KH,            TC,            QR
    A = np.array([
        [Fraction(1),   Fraction(-1, 2), Fraction(-1, 2), Fraction(0),   Fraction(0),   Fraction(0),   Fraction(0),   Fraction(0)],   # ME Eq.
        [Fraction(-1, 4), Fraction(1),   Fraction(0),   Fraction(-1, 4), Fraction(-1, 4), Fraction(-1, 4), Fraction(0),   Fraction(0)],   # GH Eq.
        [Fraction(-1, 2), Fraction(0),   Fraction(1),   Fraction(0),   Fraction(0),   Fraction(-1, 2), Fraction(0),   Fraction(0)],   # CK Eq.
        [Fraction(0),   Fraction(-1, 4), Fraction(0),   Fraction(1),   Fraction(0),   Fraction(0),   Fraction(-1, 4), Fraction(0)],   # SP Eq.
        [Fraction(0),   Fraction(-1, 3), Fraction(0),   Fraction(0),   Fraction(1),   Fraction(-1, 3), Fraction(0),   Fraction(-1, 3)],# L Eq.
        [Fraction(0),   Fraction(-1, 4), Fraction(-1, 4), Fraction(0),   Fraction(-1, 4), Fraction(1),   Fraction(0),   Fraction(-1, 4)],# KH Eq.
        [Fraction(0),   Fraction(0),   Fraction(0),   Fraction(-1, 2), Fraction(0),   Fraction(0),   Fraction(1),   Fraction(0)],   # TC Eq.
        [Fraction(0),   Fraction(0),   Fraction(0),   Fraction(0),   Fraction(-1, 3), Fraction(-1, 3), Fraction(0),   Fraction(1)]    # QR Eq.
    ], dtype=object)

    b = np.array([
        Fraction(0),      # ME: No direct path to Treasure Room
        Fraction(0),      # GH: No direct path to Treasure Room
        Fraction(0),      # CK: No direct path to Treasure Room
        Fraction(1, 4),   # SP: 1/4 chance to go directly to Treasure Room
        Fraction(0),      # L:  No direct path to Treasure Room
        Fraction(0),      # KH: No direct path to Treasure Room
        Fraction(0),      # TC: No direct path to Treasure Room
        Fraction(1, 3)    # QR: 1/3 chance to go directly to Treasure Room
    ], dtype=object)

    # Solve the system Ax = b for x
    probabilities = np.linalg.solve(A, b)

    prob_me = probabilities[0]

    # Print the equation that represents the probability for the Main Entrance
    print("The probability of winning from the Main Entrance (P_ME) is given by the equation:")
    print(f"P_ME = {A[0,1]*-1} * P_GH + {A[0,2]*-1} * P_CK")
    print(f"P_ME = 1/2 * {probabilities[1]} + 1/2 * {probabilities[2]}")
    final_calc_val = Fraction(1,2) * probabilities[1] + Fraction(1,2) * probabilities[2]
    print(f"P_ME = {final_calc_val.numerator}/{final_calc_val.denominator}")
    
    print("\n--------------------------------------------------------------")
    print("The exact probability of reaching the Treasure Room from the Main Entrance is:")
    print(f"{prob_me.numerator}/{prob_me.denominator}")

solve_bran_castle_challenge()