import numpy as np
from fractions import Fraction

def solve_bran_castle_challenge():
    """
    This function solves the Bran Castle Escape Challenge by setting up and solving
    a system of linear equations representing the probabilities of reaching the
    Treasure Room from each room in the castle.
    """

    # The system of equations is written in the form Ax = b, where x is the vector of probabilities.
    # x = [p_MainEntrance, p_GreatHall, p_CastleKitchen, p_SecretPassage, 
    #      p_Library, p_KnightsHall, p_TortureChamber, p_QueenRoom]
    
    # The coefficient matrix A is derived by rearranging the equations:
    # e.g., p_ME = 0.5*p_GH + 0.5*p_CK  =>  1*p_ME - 0.5*p_GH - 0.5*p_CK = 0
    A_frac = [
        # ME,   GH,      CK,      SP,      L,       KH,      TC,      QR
        [Fraction(1), Fraction(-1, 2), Fraction(-1, 2), Fraction(0), Fraction(0),   Fraction(0),   Fraction(0),   Fraction(0)],      # Main Entrance
        [Fraction(-1, 4), Fraction(1),   Fraction(0),   Fraction(-1, 4), Fraction(-1, 4), Fraction(-1, 4), Fraction(0),   Fraction(0)],      # Great Hall
        [Fraction(-1, 2), Fraction(0),   Fraction(1),   Fraction(0),   Fraction(0),   Fraction(-1, 2), Fraction(0),   Fraction(0)],      # Castle Kitchen
        [Fraction(0),   Fraction(-1, 4), Fraction(0),   Fraction(1),   Fraction(0),   Fraction(0),   Fraction(-1, 4), Fraction(0)],      # Secret Passage
        [Fraction(0),   Fraction(-1, 3), Fraction(0),   Fraction(0),   Fraction(1),   Fraction(-1, 3), Fraction(0),   Fraction(-1, 3)], # Library
        [Fraction(0),   Fraction(-1, 4), Fraction(-1, 4), Fraction(0),   Fraction(-1, 4), Fraction(1),   Fraction(0),   Fraction(-1, 4)], # Knights' Hall
        [Fraction(0),   Fraction(0),   Fraction(0),   Fraction(-1, 2), Fraction(0),   Fraction(0),   Fraction(1),   Fraction(0)],      # Torture Chamber
        [Fraction(0),   Fraction(0),   Fraction(0),   Fraction(0),   Fraction(-1, 3), Fraction(-1, 3), Fraction(0),   Fraction(1)]       # Queen's Room
    ]

    # The constant vector b contains the constants from the equations, which arise from
    # direct paths to the Treasure Room (probability 1) or Vampire's Lair (probability 0).
    b_frac = [
        Fraction(0),      # Main Entrance
        Fraction(0),      # Great Hall
        Fraction(0),      # Castle Kitchen
        Fraction(1, 4),   # Secret Passage (1/4 chance to go to Treasure Room)
        Fraction(0),      # Library
        Fraction(0),      # Knights' Hall
        Fraction(0),      # Torture Chamber
        Fraction(1, 3)    # Queen's Room (1/3 chance to go to Treasure Room)
    ]

    # Convert fractions to floats for NumPy's solver
    A = np.array([[float(f) for f in row] for row in A_frac])
    b = np.array([float(f) for f in b_frac])

    # Solve the system Ax = b
    solution = np.linalg.solve(A, b)

    # The probability for the Main Entrance is the first element of the solution vector
    prob_main_entrance_float = solution[0]

    # Convert the float result to an exact irreducible fraction
    prob_main_entrance_frac = Fraction(prob_main_entrance_float).limit_denominator()
    
    numerator = prob_main_entrance_frac.numerator
    denominator = prob_main_entrance_frac.denominator

    print(f"The exact probability of reaching the Treasure Room from the Main Entrance is an irreducible fraction.")
    print(f"Probability = {numerator} / {denominator}")

solve_bran_castle_challenge()
<<<109/161>>>