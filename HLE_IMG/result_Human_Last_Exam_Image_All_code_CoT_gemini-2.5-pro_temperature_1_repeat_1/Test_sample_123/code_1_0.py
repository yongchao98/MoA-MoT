import numpy as np
from fractions import Fraction

def solve_castle_escape():
    """
    Solves the Bran Castle Escape Challenge by setting up and solving a system of linear equations.
    """
    # The problem can be modeled as an absorbing Markov chain. Let p_i be the probability
    # of reaching the Treasure Room from state i. We are looking for p_ME.
    # This leads to a system of linear equations (I - Q)p = R, where Q holds the
    # transition probabilities between transient states and R holds the probabilities
    # of moving from a transient state to the treasure room.

    # State mapping to indices: 0:ME, 1:GH, 2:CK, 3:SP, 4:L, 5:KH, 6:TC, 7:QR
    
    # Matrix A = (I - Q)
    A = np.array([
        [1,    -1/2,  -1/2,  0,     0,     0,     0,     0   ],  # ME -> GH, CK
        [-1/4, 1,     0,     -1/4,  -1/4,  -1/4,  0,     0   ],  # GH -> ME, SP, L, KH
        [-1/2, 0,     1,     0,     0,     -1/2,  0,     0   ],  # CK -> ME, KH
        [0,    -1/4,  0,     1,     0,     0,     -1/4,  0   ],  # SP -> GH, TC, TR(b), VL
        [0,    -1/3,  0,     0,     1,     -1/3,  0,     -1/3],  # L  -> GH, QR, KH
        [0,    -1/4,  -1/4,  0,     -1/4,  1,     0,     -1/4],  # KH -> GH, L, QR, CK
        [0,    0,     0,     -1/2,  0,     0,     1,     0   ],  # TC -> SP, VL
        [0,    0,     0,     0,     -1/3,  -1/3,  0,     1   ]   # QR -> L, KH, TR(b)
    ])

    # Vector b = R_TR (direct probabilities to Treasure Room)
    b = np.array([0, 0, 0, 1/4, 0, 0, 0, 1/3])

    # Solve the system Ax = b for x (the vector of probabilities)
    probabilities = np.linalg.solve(A, b)

    # Extract the probabilities for the relevant rooms and convert to fractions
    p_me = Fraction(probabilities[0]).limit_denominator()
    p_gh = Fraction(probabilities[1]).limit_denominator()
    p_ck = Fraction(probabilities[2]).limit_denominator()
    
    # The probability from the Main Entrance is given by the equation:
    # P_ME = 1/2 * P_GH + 1/2 * P_CK
    # We substitute the solved values to get the final answer.
    
    print("The probability of reaching the Treasure Room from the Main Entrance, P(ME), is calculated as:")
    print("P(ME) = 1/2 * P(Great Hall) + 1/2 * P(Castle Kitchen)")
    print("\nAfter solving the system of equations for all rooms, we find the probabilities:")
    print(f"P(Great Hall) = {p_gh.numerator}/{p_gh.denominator}")
    print(f"P(Castle Kitchen) = {p_ck.numerator}/{p_ck.denominator}")
    
    print("\nSubstituting these values into the equation for P(ME):")
    print(f"P(ME) = 1/2 * {p_gh.numerator}/{p_gh.denominator} + 1/2 * {p_ck.numerator}/{p_ck.denominator}")
    
    # Perform the final calculation with fractions for exactness
    final_prob = Fraction(1, 2) * p_gh + Fraction(1, 2) * p_ck
    
    print(f"P(ME) = {final_prob.numerator}/{final_prob.denominator}")

solve_castle_escape()