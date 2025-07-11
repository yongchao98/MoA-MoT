import sympy
from sympy import symbols, Eq, solve, S

def solve_castle_escape():
    """
    This function sets up and solves the Bran Castle escape problem as a system of linear equations.
    """
    # Step 1: Define symbolic variables for the unknown probabilities for each room.
    # P_X represents the probability of eventually reaching the Treasure Room from room X.
    (P_ME, P_GH, P_CK, P_SP, 
     P_L, P_TC, P_QR, P_KH) = symbols('P_ME P_GH P_CK P_SP P_L P_TC P_QR P_KH')

    # Step 2: Define known probabilities for the absorbing states.
    # We use sympy.S() to ensure these are treated as exact rational numbers.
    P_TR = S(1)  # Probability from the Treasure Room is 1 (Success)
    P_VL = S(0)  # Probability from the Vampire's Lair is 0 (Failure)

    # Step 3: Formulate the system of 8 linear equations based on the room connections.
    
    # Eq 1: Main Entrance -> Great Hall (1/2), Castle Kitchen (1/2)
    eq1 = Eq(P_ME, S(1)/2 * P_GH + S(1)/2 * P_CK)

    # Eq 2: Great Hall -> Main Entrance (1/4), Secret Passage (1/4), Library (1/4), Knights' Hall (1/4)
    eq2 = Eq(P_GH, S(1)/4 * P_ME + S(1)/4 * P_SP + S(1)/4 * P_L + S(1)/4 * P_KH)

    # Eq 3: Castle Kitchen -> Main Entrance (1/2), Knights' Hall (1/2)
    eq3 = Eq(P_CK, S(1)/2 * P_ME + S(1)/2 * P_KH)

    # Eq 4: Secret Passage -> Great Hall (1/4), Torture Chamber (1/4), Treasure Room (1/4), Vampire's Lair (1/4)
    eq4 = Eq(P_SP, S(1)/4 * P_GH + S(1)/4 * P_TC + S(1)/4 * P_TR + S(1)/4 * P_VL)

    # Eq 5: Library -> Great Hall (1/3), Queen's Room (1/3), Knights' Hall (1/3)
    eq5 = Eq(P_L, S(1)/3 * P_GH + S(1)/3 * P_QR + S(1)/3 * P_KH)

    # Eq 6: Torture Chamber -> Secret Passage (1/2), Vampire's Lair (1/2)
    eq6 = Eq(P_TC, S(1)/2 * P_SP + S(1)/2 * P_VL)

    # Eq 7: Queen's Room -> Library (1/3), Knights' Hall (1/3), Treasure Room (1/3)
    eq7 = Eq(P_QR, S(1)/3 * P_L + S(1)/3 * P_KH + S(1)/3 * P_TR)

    # Eq 8: Knights' Hall -> Great Hall (1/4), Library (1/4), Queen's Room (1/4), Castle Kitchen (1/4)
    eq8 = Eq(P_KH, S(1)/4 * P_GH + S(1)/4 * P_L + S(1)/4 * P_QR + S(1)/4 * P_CK)
    
    print("This problem can be solved by setting up a system of linear equations.")
    print("Let P_X be the probability of reaching the Treasure Room from room X.")
    print("The system of equations is as follows (with P_TR=1 and P_VL=0 substituted):\n")
    print(f"1 * P_ME - 1/2 * P_GH - 1/2 * P_CK = 0")
    print(f"1 * P_GH - 1/4 * P_ME - 1/4 * P_SP - 1/4 * P_L - 1/4 * P_KH = 0")
    print(f"1 * P_CK - 1/2 * P_ME - 1/2 * P_KH = 0")
    print(f"1 * P_SP - 1/4 * P_GH - 1/4 * P_TC = 1/4")
    print(f"1 * P_L - 1/3 * P_GH - 1/3 * P_QR - 1/3 * P_KH = 0")
    print(f"1 * P_TC - 1/2 * P_SP = 0")
    print(f"1 * P_QR - 1/3 * P_L - 1/3 * P_KH = 1/3")
    print(f"1 * P_KH - 1/4 * P_GH - 1/4 * P_L - 1/4 * P_QR - 1/4 * P_CK = 0")

    # Step 4: Solve the system for all unknown probabilities.
    solution = solve(
        [eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8],
        [P_ME, P_GH, P_CK, P_SP, P_L, P_TC, P_QR, P_KH]
    )

    # Step 5: Extract the probability for the Main Entrance and display the result.
    prob_me = solution[P_ME]
    
    print("\nSolving this system gives the exact probability for starting at the Main Entrance.")
    print("\nThe final equation for the probability is:")
    # The result from sympy is already an irreducible fraction.
    numerator = prob_me.p
    denominator = prob_me.q
    
    print(f"P(Treasure | start at Main Entrance) = {numerator} / {denominator}")

if __name__ == '__main__':
    solve_castle_escape()