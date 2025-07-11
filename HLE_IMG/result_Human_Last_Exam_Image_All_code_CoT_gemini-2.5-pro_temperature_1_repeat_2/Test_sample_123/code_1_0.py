import numpy as np
from fractions import Fraction

def solve_bran_castle_escape():
    """
    This function sets up and solves the system of linear equations for the Bran Castle escape problem.
    """
    # Let the probabilities of reaching the Treasure Room from each location be:
    # x0 = P_ME (Main Entrance)
    # x1 = P_GH (Great Hall)
    # x2 = P_CK (Castle Kitchen)
    # x3 = P_SP (Secret Passage)
    # x4 = P_L  (Library)
    # x5 = P_KH (Knights' Hall)
    # x6 = P_QR (Queen's Room)
    # x7 = P_TC (Torture Chamber)
    # The absorbing states are Treasure Room (P_TR = 1) and Vampire's Lair (P_VL = 0).

    # We set up a system of linear equations Ax = b, where x is the vector of probabilities.
    # The equations are derived from the state transition rules:
    # P_i = sum(P_j / N_i) for all j connected to i, where N_i is the number of exits from i.
    #
    # 1. P_ME = 1/2*P_GH + 1/2*P_CK
    #    => 1*x0 - 0.5*x1 - 0.5*x2 = 0
    # 2. P_GH = 1/4*P_ME + 1/4*P_SP + 1/4*P_L + 1/4*P_KH
    #    => -0.25*x0 + 1*x1 - 0.25*x3 - 0.25*x4 - 0.25*x5 = 0
    # 3. P_CK = 1/2*P_ME + 1/2*P_KH
    #    => -0.5*x0 + 1*x2 - 0.5*x5 = 0
    # 4. P_SP = 1/4*P_GH + 1/4*P_TC + 1/4*P_TR (where P_TR=1)
    #    => -0.25*x1 + 1*x3 - 0.25*x7 = 0.25
    # 5. P_L = 1/3*P_GH + 1/3*P_QR + 1/3*P_KH
    #    => -1/3*x1 + 1*x4 - 1/3*x5 - 1/3*x6 = 0
    # 6. P_KH = 1/4*P_GH + 1/4*P_L + 1/4*P_QR + 1/4*P_CK
    #    => -0.25*x1 - 0.25*x2 - 0.25*x4 + 1*x5 - 0.25*x6 = 0
    # 7. P_QR = 1/3*P_L + 1/3*P_KH + 1/3*P_TR (where P_TR=1)
    #    => -1/3*x4 - 1/3*x5 + 1*x6 = 1/3
    # 8. P_TC = 1/2*P_SP + 1/2*P_VL (where P_VL=0)
    #    => -0.5*x3 + 1*x7 = 0

    # Coefficient matrix A
    A = np.array([
        [1,    -0.5,      -0.5,       0,         0,         0,         0,         0],
        [-0.25, 1,          0,     -0.25,     -0.25,     -0.25,         0,         0],
        [-0.5,  0,          1,         0,         0,      -0.5,         0,         0],
        [0,    -0.25,      0,         1,         0,         0,         0,     -0.25],
        [0,    -1/3,       0,         0,         1,      -1/3,      -1/3,         0],
        [0,    -0.25,    -0.25,       0,     -0.25,         1,      -0.25,        0],
        [0,       0,        0,         0,      -1/3,      -1/3,         1,         0],
        [0,       0,        0,      -0.5,         0,         0,         0,         1]
    ])

    # Constant vector b
    b = np.array([0, 0, 0, 0.25, 0, 0, 1/3, 0])

    # Solve the system Ax = b for x
    solutions = np.linalg.solve(A, b)

    # Extract the probability for the Main Entrance (x0)
    prob_main_entrance = solutions[0]

    # Convert the result to an irreducible fraction
    fraction_result = Fraction(prob_main_entrance).limit_denominator()

    # Print the final answer
    print(f"The exact probability of reaching the Treasure Room from the Main Entrance is {fraction_result.numerator}/{fraction_result.denominator}.")

    # Show the final calculation for P_MainEntrance using the solved values
    # This fulfills the request to "output each number in the final equation"
    P_GH = Fraction(solutions[1]).limit_denominator()
    P_CK = Fraction(solutions[2]).limit_denominator()

    print("\n--- Calculation Breakdown ---")
    print("The probability from the Main Entrance (P_ME) is given by the equation:")
    print(f"P_ME = (1/2) * P_GreatHall + (1/2) * P_CastleKitchen")
    print("\nAfter solving the entire system, we found the probabilities for the next rooms:")
    print(f"P_GreatHall = {P_GH.numerator}/{P_GH.denominator}")
    print(f"P_CastleKitchen = {P_CK.numerator}/{P_CK.denominator}")
    print("\nPlugging these into the equation for P_ME:")
    print(f"P_ME = (1/2) * ({P_GH.numerator}/{P_GH.denominator}) + (1/2) * ({P_CK.numerator}/{P_CK.denominator}) = {fraction_result.numerator}/{fraction_result.denominator}")

if __name__ == '__main__':
    solve_bran_castle_escape()
<<<4/17>>>