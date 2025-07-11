import numpy as np
from fractions import Fraction

def solve_bran_castle_challenge():
    """
    Solves the Bran Castle Escape Challenge by setting up and solving a system of linear equations.
    """
    
    # The rooms are the states in our system. Let p_ROOM be the probability
    # of reaching the Treasure Room from that room.
    # The absorbing states are Treasure Room (p=1) and Vampire's Lair (p=0).
    # We have 8 transient states with unknown probabilities:
    # Main Entrance (ME), Great Hall (GH), Castle Kitchen (CK), Secret Passage (SP),
    # Library (L), Knights' Hall (KH), Torture Chamber (TC), Queen's Room (QR).

    # The system of linear equations is derived from the transitions:
    # p_i = sum(P(i->j) * p_j) for all neighbors j of i.
    
    print("This script solves a system of 8 linear equations to find the probability.")
    print("The equations represent the probability of success from each room:")
    print("p_ME = 1/2 * p_GH + 1/2 * p_CK")
    print("p_GH = 1/4 * p_ME + 1/4 * p_SP + 1/4 * p_L + 1/4 * p_KH")
    print("p_CK = 1/2 * p_ME + 1/2 * p_KH")
    print("p_SP = 1/4 * p_GH + 1/4 * p_TC + 1/4 * p_TreasureRoom (1) + 1/4 * p_VampireLair (0)")
    print("p_L  = 1/3 * p_GH + 1/3 * p_KH + 1/3 * p_QR")
    print("p_KH = 1/4 * p_GH + 1/4 * p_CK + 1/4 * p_L + 1/4 * p_QR")
    print("p_TC = 1/2 * p_SP + 1/2 * p_VampireLair (0)")
    print("p_QR = 1/3 * p_L + 1/3 * p_KH + 1/3 * p_TreasureRoom (1)")
    print("-" * 30)

    # Let's map the states to indices for our matrix:
    # 0: ME, 1: GH, 2: CK, 3: SP, 4: L, 5: KH, 6: TC, 7: QR
    # We rearrange the equations into the standard form A*x = b, where x is the vector of probabilities.
    # x = [p_ME, p_GH, p_CK, p_SP, p_L, p_KH, p_TC, p_QR]^T

    # The matrix A, representing the coefficients of the variables (p_i).
    A = np.array([
        #       ME,   GH,   CK,   SP,    L,   KH,   TC,   QR
        [Fraction(1), Fraction(-1, 2), Fraction(-1, 2), Fraction(0), Fraction(0), Fraction(0), Fraction(0), Fraction(0)],       # ME Eq
        [Fraction(-1, 4), Fraction(1), Fraction(0), Fraction(-1, 4), Fraction(-1, 4), Fraction(-1, 4), Fraction(0), Fraction(0)],   # GH Eq
        [Fraction(-1, 2), Fraction(0), Fraction(1), Fraction(0), Fraction(0), Fraction(-1, 2), Fraction(0), Fraction(0)],       # CK Eq
        [Fraction(0), Fraction(-1, 4), Fraction(0), Fraction(1), Fraction(0), Fraction(0), Fraction(-1, 4), Fraction(0)],        # SP Eq
        [Fraction(0), Fraction(-1, 3), Fraction(0), Fraction(0), Fraction(1), Fraction(-1, 3), Fraction(0), Fraction(-1, 3)],   # L Eq
        [Fraction(0), Fraction(-1, 4), Fraction(-1, 4), Fraction(0), Fraction(-1, 4), Fraction(1), Fraction(0), Fraction(-1, 4)],   # KH Eq
        [Fraction(0), Fraction(0), Fraction(0), Fraction(-1, 2), Fraction(0), Fraction(0), Fraction(1), Fraction(0)],           # TC Eq
        [Fraction(0), Fraction(0), Fraction(0), Fraction(0), Fraction(-1, 3), Fraction(-1, 3), Fraction(0), Fraction(1)]        # QR Eq
    ], dtype=object)

    # The vector b, representing the constant terms from transitions to absorbing states.
    b = np.array([
        Fraction(0),      # ME
        Fraction(0),      # GH
        Fraction(0),      # CK
        Fraction(1, 4),   # SP -> Treasure Room (prob=1/4 * 1)
        Fraction(0),      # L
        Fraction(0),      # KH
        Fraction(0),      # TC -> Vampire's Lair (prob=1/2 * 0)
        Fraction(1, 3)    # QR -> Treasure Room (prob=1/3 * 1)
    ], dtype=object)

    try:
        # Solve the system of linear equations A*x = b
        solution = np.linalg.solve(A, b)
        
        # The probability from the Main Entrance is the first element of the solution vector.
        prob_me = solution[0]
        
        numerator = prob_me.numerator
        denominator = prob_me.denominator

        print("The exact probability of reaching the Treasure Room from the Main Entrance is given by the final equation:")
        print(f"p_ME = {numerator} / {denominator}")

    except np.linalg.LinAlgError:
        print("The system of equations could not be solved (the matrix is singular).")

solve_bran_castle_challenge()
<<<109/161>>>