import numpy as np
from fractions import Fraction

def solve_castle_escape():
    """
    Sets up and solves the system of linear equations for the Bran Castle escape problem.
    """
    # The system of equations can be written as Ax = b.
    # The variables are the probabilities for the 8 transient rooms:
    # 0: p_ME, 1: p_GH, 2: p_CK, 3: p_SP, 4: p_L, 5: p_KH, 6: p_TC, 7: p_QR

    # A is the matrix of coefficients from the rearranged equations (e.g., p_ME - 0.5*p_GH - 0.5*p_CK = 0)
    A = np.array([
        [1,      -1/2,   -1/2,   0,      0,      0,      0,      0   ],   # Main Entrance
        [-1/4,   1,      0,      -1/4,   -1/4,   -1/4,   0,      0   ],   # Great Hall
        [-1/2,   0,      1,      0,      0,      -1/2,   0,      0   ],   # Castle Kitchen
        [0,      -1/4,   0,      1,      0,      0,      -1/4,   0   ],   # Secret Passage
        [0,      -1/3,   0,      0,      1,      -1/3,   0,      -1/3],   # Library
        [0,      -1/4,   -1/4,   0,      -1/4,   1,      0,      -1/4],   # Knights' Hall
        [0,      0,      0,      -1/2,   0,      0,      1,      0   ],   # Torture Chamber
        [0,      0,      0,      0,      -1/3,   -1/3,   0,      1   ]    # Queen's Room
    ])

    # b is the vector of constants, which come from direct paths to the Treasure Room.
    b = np.array([
        0,          # Main Entrance
        0,          # Great Hall
        0,          # Castle Kitchen
        1/4,        # Secret Passage
        0,          # Library
        0,          # Knights' Hall
        0,          # Torture Chamber
        1/3         # Queen's Room
    ])

    # Solve the system of linear equations
    try:
        probabilities = np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        print("The system of equations could not be solved.")
        return

    # The probability of reaching the Treasure Room from the Main Entrance is the first element
    p_main_entrance = probabilities[0]

    # Convert the floating-point result to an exact irreducible fraction
    final_fraction = Fraction(p_main_entrance).limit_denominator()
    
    numerator = final_fraction.numerator
    denominator = final_fraction.denominator

    print(f"The probability of eventually reaching the Treasure Room from the Main Entrance is:")
    print(f"The final equation is the irreducible fraction: {numerator} / {denominator}")

solve_castle_escape()