from sympy import symbols, Eq, solve, S

def solve_bran_castle_challenge():
    """
    This function sets up and solves the system of linear equations
    for the Bran Castle escape problem to find the exact probability
    of reaching the Treasure Room from the Main Entrance.
    """
    # Define the probabilities for each room as symbolic variables.
    # P_XX represents the probability of reaching the treasure room from room XX.
    P_ME, P_GH, P_CK, P_SP, P_L, P_KH, P_TC, P_QR = symbols('P_ME P_GH P_CK P_SP P_L P_KH P_TC P_QR')

    # The probabilities for the absorbing states are known.
    # S() converts numbers to sympy's rational type for exact calculations.
    P_TR = S(1)  # Probability of success if starting in the Treasure Room is 1.
    P_VL = S(0)  # Probability of success if starting in the Vampire's Lair is 0.

    # Based on the castle layout, we can set up one equation per room.
    # The probability from a room is the weighted average of the probabilities of the rooms it leads to.

    # 1. Main Entrance (ME) -> Great Hall (GH), Castle Kitchen (CK)
    eq_ME = Eq(P_ME, S('1/2') * P_GH + S('1/2') * P_CK)

    # 2. Great Hall (GH) -> Main Entrance (ME), Secret Passage (SP), Library (L), Knights' Hall (KH)
    eq_GH = Eq(P_GH, S('1/4') * P_ME + S('1/4') * P_SP + S('1/4') * P_L + S('1/4') * P_KH)

    # 3. Castle Kitchen (CK) -> Main Entrance (ME), Knights' Hall (KH)
    eq_CK = Eq(P_CK, S('1/2') * P_ME + S('1/2') * P_KH)

    # 4. Secret Passage (SP) -> Great Hall (GH), Torture Chamber (TC), Treasure Room (TR), Vampire's Lair (VL)
    eq_SP = Eq(P_SP, S('1/4') * P_GH + S('1/4') * P_TC + S('1/4') * P_TR + S('1/4') * P_VL)

    # 5. Library (L) -> Great Hall (GH), Queen's Room (QR), Knights' Hall (KH)
    eq_L = Eq(P_L, S('1/3') * P_GH + S('1/3') * P_QR + S('1/3') * P_KH)

    # 6. Knights' Hall (KH) -> Great Hall (GH), Library (L), Queen's Room (QR), Castle Kitchen (CK)
    eq_KH = Eq(P_KH, S('1/4') * P_GH + S('1/4') * P_L + S('1/4') * P_QR + S('1/4') * P_CK)

    # 7. Torture Chamber (TC) -> Secret Passage (SP), Vampire's Lair (VL)
    eq_TC = Eq(P_TC, S('1/2') * P_SP + S('1/2') * P_VL)

    # 8. Queen's Room (QR) -> Library (L), Knights' Hall (KH), Treasure Room (TR)
    eq_QR = Eq(P_QR, S('1/3') * P_L + S('1/3') * P_KH + S('1/3') * P_TR)

    # Solve the system of 8 linear equations for the 8 unknown probabilities
    solution = solve(
        (eq_ME, eq_GH, eq_CK, eq_SP, eq_L, eq_KH, eq_TC, eq_QR),
        (P_ME, P_GH, P_CK, P_SP, P_L, P_KH, P_TC, P_QR)
    )

    # The result we are interested in is the probability from the Main Entrance (P_ME)
    prob_ME = solution[P_ME]
    
    # Extract the numerator and denominator to print the final equation
    numerator = prob_ME.p
    denominator = prob_ME.q
    
    print(f"The exact probability of reaching the Treasure Room from the Main Entrance is given by the equation:")
    print(f"P(Success | Start at Main Entrance) = {numerator} / {denominator}")

if __name__ == '__main__':
    solve_bran_castle_challenge()