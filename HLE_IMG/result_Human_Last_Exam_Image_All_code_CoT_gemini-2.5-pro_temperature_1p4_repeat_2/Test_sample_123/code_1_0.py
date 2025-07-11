import numpy as np
from fractions import Fraction

def solve_bran_castle_escape():
    """
    Calculates the probability of reaching the Treasure Room from the Main Entrance.

    This function models the castle as a system of linear equations where each
    variable represents the probability of reaching the treasure from a specific room.
    It uses numpy to solve this system with exact fractional arithmetic.
    """
    
    # Let the variables for the probabilities of reaching the treasure from each room be ordered as:
    # 0: Main Entrance (me)
    # 1: Great Hall (gh)
    # 2: Castle Kitchen (ck)
    # 3: Secret Passage (sp)
    # 4: Library (l)
    # 5: Torture Chamber (tc)
    # 6: Queen's Room (qr)
    # 7: Knights' Hall (kh)

    # We set up the system of equations in the form Ax = b.
    # The matrix A contains the coefficients of the probability variables.
    # The vector b contains the constant terms from the equations.

    # A*x = b, where x = [p_me, p_gh, p_ck, p_sp, p_l, p_tc, p_qr, p_kh]^T
    # The equations are rearranged to the form: c1*p_me + c2*p_gh + ... = k

    A = np.array([
        #      me,       gh,       ck,       sp,       l,        tc,       qr,       kh
        [Fraction(1), Fraction(-1, 2), Fraction(-1, 2), Fraction(0),   Fraction(0),    Fraction(0),    Fraction(0),    Fraction(0)],    # eq_me: p_me - 1/2*p_gh - 1/2*p_ck = 0
        [Fraction(-1, 4), Fraction(1),    Fraction(0),    Fraction(-1, 4), Fraction(-1, 4), Fraction(0),    Fraction(0),    Fraction(-1, 4)],  # eq_gh
        [Fraction(-1, 2), Fraction(0),    Fraction(1),    Fraction(0),    Fraction(0),    Fraction(0),    Fraction(0),    Fraction(-1, 2)],  # eq_ck
        [Fraction(0),   Fraction(-1, 4), Fraction(0),    Fraction(1),    Fraction(0),    Fraction(-1, 4), Fraction(0),    Fraction(0)],    # eq_sp
        [Fraction(0),   Fraction(-1, 3), Fraction(0),    Fraction(0),    Fraction(1),    Fraction(0),    Fraction(-1, 3), Fraction(-1, 3)],  # eq_l
        [Fraction(0),   Fraction(0),    Fraction(0),    Fraction(-1, 2), Fraction(0),    Fraction(1),    Fraction(0),    Fraction(0)],    # eq_tc
        [Fraction(0),   Fraction(0),    Fraction(0),    Fraction(0),    Fraction(-1, 3), Fraction(0),    Fraction(1),    Fraction(-1, 3)],  # eq_qr
        [Fraction(0),   Fraction(-1, 4), Fraction(-1, 4), Fraction(0),    Fraction(-1, 4), Fraction(0),    Fraction(-1, 4), Fraction(1)],    # eq_kh
    ], dtype=object)

    b = np.array([
        Fraction(0),      # eq_me constant
        Fraction(0),      # eq_gh constant
        Fraction(0),      # eq_ck constant
        Fraction(1, 4),   # eq_sp constant (from 1/4 * P(Treasure Room))
        Fraction(0),      # eq_l constant
        Fraction(0),      # eq_tc constant
        Fraction(1, 3),   # eq_qr constant (from 1/3 * P(Treasure Room))
        Fraction(0)       # eq_kh constant
    ], dtype=object)

    # Solve the system of linear equations
    solution = np.linalg.solve(A, b)

    # The probability of reaching the treasure from the Main Entrance is the first element
    p_main_entrance = solution[0]
    
    # Print the final equation P(Main Entrance) = numerator / denominator
    print("The final probability equation is:")
    print(f"P(Main Entrance) = {p_main_entrance.numerator} / {p_main_entrance.denominator}")

solve_bran_castle_escape()