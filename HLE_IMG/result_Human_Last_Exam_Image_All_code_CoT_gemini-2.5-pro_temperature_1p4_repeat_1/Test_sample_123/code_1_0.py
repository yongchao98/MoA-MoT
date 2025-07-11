import numpy as np
from fractions import Fraction

def solve_bran_castle_challenge():
    """
    Solves the Bran Castle escape challenge by setting up and solving a system of
    linear equations representing the probabilities of reaching the Treasure Room
    from each room in the castle.
    """
    
    # The transient states (non-terminal rooms) are ordered as follows:
    # 0: Main Entrance (ME)
    # 1: Great Hall (GH)
    # 2: Castle Kitchen (CK)
    # 3: Secret Passage (SP)
    # 4: Library (L)
    # 5: Knights' Hall (KH)
    # 6: Queen's Room (QR)
    # 7: Torture Chamber (TC)
    
    # We want to solve a system of linear equations of the form Ax = b,
    # where x is the vector of probabilities for each room.
    # The matrix A is (I - T), where I is the identity matrix and T is the
    # transition matrix between transient states. The vector b represents the
    # probability of moving directly to the Treasure Room (success state).

    # Equations derived from the problem description:
    # P_ME = (1/2)*P_GH + (1/2)*P_CK
    # P_GH = (1/4)*P_ME + (1/4)*P_SP + (1/4)*P_L + (1/4)*P_KH
    # P_CK = (1/2)*P_ME + (1/2)*P_KH
    # P_SP = (1/4)*P_GH + (1/4)*P_TC + (1/4)*P_TR + (1/4)*P_VL  => P_SP - (1/4)P_GH - (1/4)P_TC = 1/4
    # P_L  = (1/3)*P_GH + (1/3)*P_QR + (1/3)*P_KH
    # P_KH = (1/4)*P_GH + (1/4)*P_L + (1/4)*P_QR + (1/4)*P_CK
    # P_QR = (1/3)*P_L + (1/3)*P_KH + (1/3)*P_TR               => P_QR - (1/3)P_L - (1/3)P_KH = 1/3
    # P_TC = (1/2)*P_SP + (1/2)*P_VL                          => P_TC - (1/2)P_SP = 0

    # A matrix for the system (I-T)x = b
    A = np.array([
        # ME,   GH,   CK,   SP,    L,   KH,   QR,   TC  <- Variables
        [ 1,  -1/2, -1/2,    0,    0,    0,    0,    0],  # ME equation
        [-1/4,  1,    0,  -1/4, -1/4, -1/4,    0,    0],  # GH equation
        [-1/2,  0,    1,    0,    0,  -1/2,    0,    0],  # CK equation
        [ 0,  -1/4,   0,    1,    0,    0,    0,  -1/4],  # SP equation
        [ 0,  -1/3,   0,    0,    1,  -1/3, -1/3,    0],  # L equation
        [ 0,  -1/4, -1/4,   0,  -1/4,   1,  -1/4,    0],  # KH equation
        [ 0,    0,    0,    0,  -1/3, -1/3,   1,     0],  # QR equation
        [ 0,    0,    0,  -1/2,   0,    0,    0,     1]   # TC equation
    ])

    # b vector (probabilities of one-step transition to Treasure Room)
    b = np.array([0, 0, 0, 1/4, 0, 0, 1/3, 0])

    # Solve the system Ax = b for x
    try:
        probabilities = np.linalg.solve(A, b)
        
        # The probability for the Main Entrance is the first element
        prob_me = probabilities[0]
        
        # Convert the floating point result to an exact fraction
        final_fraction = Fraction(prob_me).limit_denominator()
        
        print("The probability of reaching the Treasure Room from the Main Entrance has been calculated.")
        print("The final equation for this probability is:")
        
        # Output each number in the final equation as requested
        numerator = final_fraction.numerator
        denominator = final_fraction.denominator
        
        print(f"P(Success) = {numerator} / {denominator}")

    except np.linalg.LinAlgError:
        print("The system of equations could not be solved. The matrix may be singular.")

solve_bran_castle_challenge()