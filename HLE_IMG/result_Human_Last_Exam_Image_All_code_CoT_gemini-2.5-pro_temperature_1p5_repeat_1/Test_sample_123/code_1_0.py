import numpy as np
from fractions import Fraction

def solve_castle_escape():
    """
    Solves the Bran Castle escape probability problem.
    """
    # Assign an index to each transient state (room)
    rooms = {
        'ME': 0, 'GH': 1, 'CK': 2, 'SP': 3,
        'L': 4, 'KH': 5, 'QR': 6, 'TC': 7
    }
    num_rooms = len(rooms)

    # Initialize the matrix A and vector b for the system A*x = b
    # where x is the vector of probabilities for each room.
    # We use Fraction for exact arithmetic.
    A = [[Fraction(0) for _ in range(num_rooms)] for _ in range(num_rooms)]
    b = [Fraction(0) for _ in range(num_rooms)]

    # Populate the A matrix and b vector based on the problem description.
    # Each row corresponds to an equation P_room = sum(prob_i * P_next_room_i)
    # which we rearrange to P_room - sum(...) = direct_prob_to_treasure.

    # Equation for Main Entrance (ME) -> GH, CK
    me_idx = rooms['ME']
    A[me_idx][rooms['ME']] = Fraction(1)
    A[me_idx][rooms['GH']] = -Fraction(1, 2)
    A[me_idx][rooms['CK']] = -Fraction(1, 2)
    b[me_idx] = Fraction(0)

    # Equation for Great Hall (GH) -> ME, SP, L, KH
    gh_idx = rooms['GH']
    A[gh_idx][rooms['GH']] = Fraction(1)
    A[gh_idx][rooms['ME']] = -Fraction(1, 4)
    A[gh_idx][rooms['SP']] = -Fraction(1, 4)
    A[gh_idx][rooms['L']] = -Fraction(1, 4)
    A[gh_idx][rooms['KH']] = -Fraction(1, 4)
    b[gh_idx] = Fraction(0)
    
    # Equation for Castle Kitchen (CK) -> ME, KH
    ck_idx = rooms['CK']
    A[ck_idx][rooms['CK']] = Fraction(1)
    A[ck_idx][rooms['ME']] = -Fraction(1, 2)
    A[ck_idx][rooms['KH']] = -Fraction(1, 2)
    b[ck_idx] = Fraction(0)

    # Equation for Secret Passage (SP) -> GH, TC, Treasure Room, Vampire's Lair
    sp_idx = rooms['SP']
    A[sp_idx][rooms['SP']] = Fraction(1)
    A[sp_idx][rooms['GH']] = -Fraction(1, 4)
    A[sp_idx][rooms['TC']] = -Fraction(1, 4)
    b[sp_idx] = Fraction(1, 4)  # 1/4 chance to go directly to Treasure Room

    # Equation for Library (L) -> GH, QR, KH
    l_idx = rooms['L']
    A[l_idx][rooms['L']] = Fraction(1)
    A[l_idx][rooms['GH']] = -Fraction(1, 3)
    A[l_idx][rooms['QR']] = -Fraction(1, 3)
    A[l_idx][rooms['KH']] = -Fraction(1, 3)
    b[l_idx] = Fraction(0)

    # Equation for Knights' Hall (KH) -> GH, L, QR, CK
    kh_idx = rooms['KH']
    A[kh_idx][rooms['KH']] = Fraction(1)
    A[kh_idx][rooms['GH']] = -Fraction(1, 4)
    A[kh_idx][rooms['L']] = -Fraction(1, 4)
    A[kh_idx][rooms['QR']] = -Fraction(1, 4)
    A[kh_idx][rooms['CK']] = -Fraction(1, 4)
    b[kh_idx] = Fraction(0)

    # Equation for Queen's Room (QR) -> L, KH, Treasure Room
    qr_idx = rooms['QR']
    A[qr_idx][rooms['QR']] = Fraction(1)
    A[qr_idx][rooms['L']] = -Fraction(1, 3)
    A[qr_idx][rooms['KH']] = -Fraction(1, 3)
    b[qr_idx] = Fraction(1, 3)  # 1/3 chance to go directly to Treasure Room

    # Equation for Torture Chamber (TC) -> SP, Vampire's Lair
    tc_idx = rooms['TC']
    A[tc_idx][rooms['TC']] = Fraction(1)
    A[tc_idx][rooms['SP']] = -Fraction(1, 2)
    b[tc_idx] = Fraction(0)
    
    # Convert Fractions to floats for numpy solver
    A_float = np.array(A, dtype=float)
    b_float = np.array(b, dtype=float)
    
    # Solve the system of linear equations
    solution_float = np.linalg.solve(A_float, b_float)
    
    # Convert float solutions back to Fractions to get the exact answer
    solution = [Fraction(x).limit_denominator() for x in solution_float]
    
    p_me = solution[rooms['ME']]
    p_gh = solution[rooms['GH']]
    p_ck = solution[rooms['CK']]

    # The fundamental equation for the starting room (Main Entrance) is:
    # P_ME = (1/2) * P_GH + (1/2) * P_CK
    # We print the calculation using the solved values.
    
    c1 = Fraction(1, 2)
    c2 = Fraction(1, 2)
    
    print("The probability of reaching the Treasure Room from the Main Entrance is calculated as:")
    print("P_ME = 1/2 * P_GH + 1/2 * P_CK")
    print("\nAfter solving the system of equations for all rooms, we found:")
    print(f"P_GH = {p_gh.numerator}/{p_gh.denominator}")
    print(f"P_CK = {p_ck.numerator}/{p_ck.denominator}")
    
    print("\nSubstituting these values into the equation for P_ME:")
    final_eq = (f"{c1.numerator}/{c1.denominator} * {p_gh.numerator}/{p_gh.denominator} + "
                f"{c2.numerator}/{c2.denominator} * {p_ck.numerator}/{p_ck.denominator} = "
                f"{p_me.numerator}/{p_me.denominator}")

    print(final_eq)

solve_castle_escape()
<<<109/161>>>