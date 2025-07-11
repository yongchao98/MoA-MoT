import numpy as np
from fractions import Fraction

def solve_bran_castle_escape():
    """
    This function solves the Bran Castle escape problem by setting up and solving
    a system of linear equations representing the probabilities of reaching the
    Treasure Room from each location.
    """

    # Let the probabilities of reaching the Treasure Room from each room be:
    # p_ME: Main Entrance (index 0)
    # p_GH: Great Hall (index 1)
    # p_CK: Castle Kitchen (index 2)
    # p_SP: Secret Passage (index 3)
    # p_L:  Library (index 4)
    # p_KH: Knights' Hall (index 5)
    # p_TC: Torture Chamber (index 6)
    # p_QR: Queen's Room (index 7)

    # The probability of reaching the Treasure Room from the Treasure Room is 1.
    # The probability of reaching the Treasure Room from the Vampire's Lair is 0.

    # We set up a system of linear equations Ax = b.
    # For any room X, p_X = sum(prob(X -> Y) * p_Y) for all adjacent rooms Y.
    # We rearrange this to the form ... = constant.

    # The equations are:
    # 1. p_ME = (1/2)p_GH + (1/2)p_CK
    #    => 1*p_ME - (1/2)p_GH - (1/2)p_CK = 0
    # 2. p_GH = (1/4)p_ME + (1/4)p_SP + (1/4)p_L + (1/4)p_KH
    #    => -(1/4)p_ME + 1*p_GH - (1/4)p_SP - (1/4)p_L - (1/4)p_KH = 0
    # 3. p_CK = (1/2)p_ME + (1/2)p_KH
    #    => -(1/2)p_ME + 1*p_CK - (1/2)p_KH = 0
    # 4. p_SP = (1/4)p_GH + (1/4)p_TC + (1/4)*p_TR + (1/4)*p_VL
    #    p_SP = (1/4)p_GH + (1/4)p_TC + (1/4)*1 + (1/4)*0
    #    => -(1/4)p_GH + 1*p_SP - (1/4)p_TC = 1/4
    # 5. p_L  = (1/3)p_GH + (1/3)p_QR + (1/3)p_KH
    #    => -(1/3)p_GH + 1*p_L - (1/3)p_KH - (1/3)p_QR = 0
    # 6. p_KH = (1/4)p_GH + (1/4)p_L + (1/4)p_QR + (1/4)p_CK
    #    => -(1/4)p_GH - (1/4)p_CK - (1/4)p_L + 1*p_KH - (1/4)p_QR = 0
    # 7. p_TC = (1/2)p_SP + (1/2)*p_VL
    #    p_TC = (1/2)p_SP + (1/2)*0
    #    => -(1/2)p_SP + 1*p_TC = 0
    # 8. p_QR = (1/3)p_L + (1/3)p_KH + (1/3)*p_TR
    #    p_QR = (1/3)p_L + (1/3)p_KH + (1/3)*1
    #    => -(1/3)p_L - (1/3)p_KH + 1*p_QR = 1/3

    # Matrix A (coefficients of the variables)
    A = np.array([
        [1,    -1/2,  -1/2,    0,     0,     0,     0,     0   ], # Eq 1
        [-1/4,   1,     0,   -1/4,  -1/4,  -1/4,    0,     0   ], # Eq 2
        [-1/2,   0,     1,     0,     0,   -1/2,    0,     0   ], # Eq 3
        [0,    -1/4,    0,     1,     0,     0,   -1/4,    0   ], # Eq 4
        [0,    -1/3,    0,     0,     1,   -1/3,    0,   -1/3  ], # Eq 5
        [0,    -1/4,  -1/4,    0,   -1/4,    1,     0,   -1/4  ], # Eq 6
        [0,      0,     0,   -1/2,    0,     0,     1,     0   ], # Eq 7
        [0,      0,     0,     0,   -1/3,  -1/3,    0,     1   ]  # Eq 8
    ])

    # Vector b (constants on the right side of the equations)
    b = np.array([0, 0, 0, 1/4, 0, 0, 0, 1/3])

    # Solve the system Ax = b for x (the vector of probabilities)
    probabilities = np.linalg.solve(A, b)

    # The probability of success from the Main Entrance is the first element
    prob_main_entrance = probabilities[0]

    # Convert the result to an irreducible fraction
    fraction_result = Fraction(prob_main_entrance).limit_denominator()

    # Print the explanation and the final result
    print("The problem is solved by setting up a system of linear equations.")
    print("Let p_X be the probability of success starting from room X.")
    print("The equations for the probabilities are:")
    print("p_ME = (1/2) * p_GH + (1/2) * p_CK")
    print("p_GH = (1/4) * p_ME + (1/4) * p_SP + (1/4) * p_L + (1/4) * p_KH")
    print("p_CK = (1/2) * p_ME + (1/2) * p_KH")
    print("p_SP = (1/4) * p_GH + (1/4) * p_TC + (1/4) * 1 + (1/4) * 0")
    print("p_L = (1/3) * p_GH + (1/3) * p_QR + (1/3) * p_KH")
    print("p_KH = (1/4) * p_GH + (1/4) * p_L + (1/4) * p_QR + (1/4) * p_CK")
    print("p_TC = (1/2) * p_SP + (1/2) * 0")
    print("p_QR = (1/3) * p_L + (1/3) * p_KH + (1/3) * 1")
    print("\nSolving this system gives the probability for the Main Entrance.")
    print(f"\nThe exact probability of reaching the Treasure Room from the Main Entrance is: {fraction_result.numerator}/{fraction_result.denominator}")

if __name__ == "__main__":
    solve_bran_castle_escape()