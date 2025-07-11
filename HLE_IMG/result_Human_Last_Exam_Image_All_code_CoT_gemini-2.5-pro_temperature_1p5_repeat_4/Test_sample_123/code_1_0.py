import numpy as np
from fractions import Fraction

def solve_bran_castle_escape():
    """
    Calculates the probability of reaching the Treasure Room from the Main Entrance
    by modeling the castle as a system of linear equations and solving for the probabilities.
    """
    # Set up the system of linear equations Ax = b, where x is the vector of 
    # probabilities of reaching the Treasure Room from each room.
    # The order of rooms (variables) is:
    # 0: Main Entrance (ME), 1: Great Hall (GH), 2: Castle Kitchen (CK), 
    # 3: Secret Passage (SP), 4: Library (L), 5: Knights' Hall (KH), 
    # 6: Torture Chamber (TC), 7: Queen's Room (QR)

    # We know P_TreasureRoom = 1 and P_VampireLair = 0.
    # The equations below are derived from P_i = sum(P_j / N_i) over all neighbors j,
    # where N_i is the number of exits from room i.

    # P_ME = (1/2)P_GH + (1/2)P_CK                       =>  2*P_ME - P_GH - P_CK = 0
    # P_GH = (1/4)P_ME + (1/4)P_SP + (1/4)P_L + (1/4)P_KH => -P_ME + 4*P_GH - P_SP - P_L - P_KH = 0
    # P_CK = (1/2)P_ME + (1/2)P_KH                       => -P_ME + 2*P_CK - P_KH = 0
    # P_SP = (1/4)P_GH + (1/4)P_TC + 1/4                  => -P_GH + 4*P_SP - P_TC = 1
    # P_L  = (1/3)P_GH + (1/3)P_QR + (1/3)P_KH           => -P_GH + 3*P_L - P_KH - P_QR = 0
    # P_KH = (1/4)P_GH + (1/4)P_L + (1/4)P_QR + (1/4)P_CK => -P_CK - P_GH - P_L + 4*P_KH - P_QR = 0
    # P_TC = (1/2)P_SP                                   => -P_SP + 2*P_TC = 0
    # P_QR = (1/3)P_L + (1/3)P_KH + 1/3                  => -P_L - P_KH + 3*P_QR = 1

    # Coefficient matrix A
    A = np.array([
        [ 2, -1, -1,  0,  0,  0,  0,  0],  # ME equation
        [-1,  4,  0, -1, -1, -1,  0,  0],  # GH equation
        [-1,  0,  2,  0,  0, -1,  0,  0],  # CK equation
        [ 0, -1,  0,  4,  0,  0, -1,  0],  # SP equation
        [ 0, -1,  0,  0,  3, -1,  0, -1],  # L equation
        [ 0, -1, -1,  0, -1,  4,  0, -1],  # KH equation
        [ 0,  0,  0, -1,  0,  0,  2,  0],  # TC equation
        [ 0,  0,  0,  0, -1, -1,  0,  3]   # QR equation
    ])

    # Constant vector b
    b = np.array([0, 0, 0, 1, 0, 0, 0, 1])

    # Solve the system of equations
    probabilities = np.linalg.solve(A, b)

    # The probability from the Main Entrance is the first element
    prob_main_entrance = probabilities[0]

    # Convert the result to an irreducible fraction
    final_fraction = Fraction(prob_main_entrance).limit_denominator()
    numerator = final_fraction.numerator
    denominator = final_fraction.denominator

    # Print the final result and the numbers in the final fraction equation
    print(f"The probability of reaching the Treasure Room from the Main Entrance is given by the fraction P = {numerator} / {denominator}.")
    print("\nHere are the numbers in the final equation:")
    print(f"Numerator: {numerator}")
    print(f"Denominator: {denominator}")

solve_bran_castle_escape()