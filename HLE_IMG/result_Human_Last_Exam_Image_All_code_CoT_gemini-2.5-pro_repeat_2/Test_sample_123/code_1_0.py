import numpy as np
from fractions import Fraction

def solve_bran_castle_escape():
    """
    Solves the Bran Castle Escape problem by setting up and solving a system of linear equations.
    """
    # The order of rooms (variables) is:
    # 0: Main Entrance (ME), 1: Great Hall (GH), 2: Castle Kitchen (CK),
    # 3: Secret Passage (SP), 4: Library (L), 5: Knights' Hall (KH),
    # 6: Torture Chamber (TC), 7: Queen's Room (QR)

    # We set up a system of linear equations Ax = b, where x is the vector of probabilities.
    # Each row in matrix A and vector b corresponds to an equation for a room.
    # The equations are rearranged to the form: c1*P1 + c2*P2 + ... = k

    # Coefficient matrix A
    A = np.array([
        # ME,  GH,  CK,  SP,   L,  KH,  TC,  QR
        [ 1, -1/2, -1/2,    0,    0,    0,   0,    0], # P_ME - 0.5*P_GH - 0.5*P_CK = 0
        [-1/4,   1,    0, -1/4, -1/4, -1/4,   0,    0], # P_GH - 0.25*P_ME - 0.25*P_SP - ... = 0
        [-1/2,   0,    1,    0,    0, -1/2,   0,    0], # P_CK - 0.5*P_ME - 0.5*P_KH = 0
        [   0,-1/4,    0,    1,    0,    0,-1/4,    0], # P_SP - 0.25*P_GH - 0.25*P_TC = 0.25
        [   0,-1/3,    0,    0,    1, -1/3,   0, -1/3], # P_L - (1/3)*P_GH - (1/3)*P_KH - ... = 0
        [   0,-1/4, -1/4,    0, -1/4,    1,   0, -1/4], # P_KH - 0.25*P_GH - 0.25*P_CK - ... = 0
        [   0,   0,    0, -1/2,    0,    0,   1,    0], # P_TC - 0.5*P_SP = 0
        [   0,   0,    0,    0, -1/3, -1/3,   0,    1]  # P_QR - (1/3)*P_L - (1/3)*P_KH = 1/3
    ])

    # Constants vector b
    # These come from connections to the Treasure Room (prob=1) or Vampire's Lair (prob=0).
    b = np.array([
        0,       # Main Entrance equation
        0,       # Great Hall equation
        0,       # Castle Kitchen equation
        1/4,     # Secret Passage equation (from 1/4 * P_TR)
        0,       # Library equation
        0,       # Knights' Hall equation
        0,       # Torture Chamber equation
        1/3      # Queen's Room equation (from 1/3 * P_TR)
    ])

    # Solve the system of equations
    probabilities = np.linalg.solve(A, b)

    # Extract the probabilities for the relevant rooms
    p_me = probabilities[0]
    p_gh = probabilities[1]
    p_ck = probabilities[2]

    # Convert the final result to an irreducible fraction
    fraction_result = Fraction(p_me).limit_denominator()
    numerator = fraction_result.numerator
    denominator = fraction_result.denominator

    print("The probability of reaching the Treasure Room from the Main Entrance is calculated based on the probabilities of the connecting rooms.")
    print("The governing equation is:")
    print(f"P(Main Entrance) = (1/2) * P(Great Hall) + (1/2) * P(Castle Kitchen)")
    print("\nAfter solving the full system of equations, we find the probabilities for the Great Hall and Castle Kitchen:")
    print(f"P(Great Hall) = {p_gh:.4f}")
    print(f"P(Castle Kitchen) = {p_ck:.4f}")
    print("\nPlugging these values into the equation for the Main Entrance:")
    print(f"P(Main Entrance) = (1 / 2) * {p_gh:.4f} + (1 / 2) * {p_ck:.4f} = {p_me:.4f}")

    print("\nTherefore, the exact probability of reaching the Treasure Room starting from the Main Entrance is:")
    print(f"{numerator} / {denominator}")

solve_bran_castle_escape()
<<<3/13>>>