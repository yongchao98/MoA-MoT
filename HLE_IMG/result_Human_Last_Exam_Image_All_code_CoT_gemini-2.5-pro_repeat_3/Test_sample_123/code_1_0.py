import sympy

def solve_bran_castle_challenge():
    """
    This function sets up and solves the system of linear equations for the Bran Castle
    Escape Challenge to find the exact probability of reaching the Treasure Room.
    """
    
    # Define the symbols for the unknown probabilities for each room.
    # P_X represents the probability of eventually reaching the Treasure Room from room X.
    P_ME, P_GH, P_CK, P_SP, P_L, P_KH, P_QR, P_TC = sympy.symbols(
        'P_ME P_GH P_CK P_SP P_L P_KH P_QR P_TC'
    )
    
    # List of variables for the solver
    variables = [P_ME, P_GH, P_CK, P_SP, P_L, P_KH, P_QR, P_TC]

    # Define the system of linear equations based on the castle layout.
    # Each equation is of the form: P_i = sum(prob(i -> j) * P_j)
    # The probability of moving to any exit is 1 / (number of exits).
    # We set P_TreasureRoom = 1 and P_VampireLair = 0.
    equations = [
        # 1. Main Entrance (2 exits) -> Great Hall, Castle Kitchen
        sympy.Eq(P_ME, sympy.S(1)/2 * P_GH + sympy.S(1)/2 * P_CK),

        # 2. Great Hall (4 exits) -> Main Entrance, Secret Passage, Library, Knights' Hall
        sympy.Eq(P_GH, sympy.S(1)/4 * P_ME + sympy.S(1)/4 * P_SP + sympy.S(1)/4 * P_L + sympy.S(1)/4 * P_KH),

        # 3. Castle Kitchen (2 exits) -> Main Entrance, Knights' Hall
        sympy.Eq(P_CK, sympy.S(1)/2 * P_ME + sympy.S(1)/2 * P_KH),

        # 4. Secret Passage (4 exits) -> Great Hall, Torture Chamber, Treasure Room (P=1), Vampire's Lair (P=0)
        sympy.Eq(P_SP, sympy.S(1)/4 * P_GH + sympy.S(1)/4 * P_TC + sympy.S(1)/4 * 1 + sympy.S(1)/4 * 0),

        # 5. Library (3 exits) -> Great Hall, Queen's Room, Knights' Hall
        sympy.Eq(P_L, sympy.S(1)/3 * P_GH + sympy.S(1)/3 * P_QR + sympy.S(1)/3 * P_KH),

        # 6. Knights' Hall (4 exits) -> Great Hall, Library, Queen's Room, Castle Kitchen
        sympy.Eq(P_KH, sympy.S(1)/4 * P_GH + sympy.S(1)/4 * P_L + sympy.S(1)/4 * P_QR + sympy.S(1)/4 * P_CK),

        # 7. Queen's Room (3 exits) -> Library, Knights' Hall, Treasure Room (P=1)
        sympy.Eq(P_QR, sympy.S(1)/3 * P_L + sympy.S(1)/3 * P_KH + sympy.S(1)/3 * 1),

        # 8. Torture Chamber (2 exits) -> Secret Passage, Vampire's Lair (P=0)
        sympy.Eq(P_TC, sympy.S(1)/2 * P_SP + sympy.S(1)/2 * 0)
    ]

    # Solve the system of equations. The result is a dictionary mapping variables to solutions.
    solution = sympy.solve(equations, variables, dict=True)
    if not solution:
        print("The system of equations could not be solved.")
        return
        
    solution = solution[0]

    # Extract the required probabilities for the final calculation display
    p_me_val = solution[P_ME]
    p_gh_val = solution[P_GH]
    p_ck_val = solution[P_CK]

    # Display the final calculation for P_ME as requested
    print("The probability of reaching the Treasure Room from the Main Entrance is determined by the probabilities of the rooms it connects to.")
    print("The governing equation is:")
    print("P_MainEntrance = (1/2) * P_GreatHall + (1/2) * P_CastleKitchen")
    print("\nSubstituting the solved probabilities:")
    print(f"P_MainEntrance = (1/2) * ({p_gh_val}) + (1/2) * ({p_ck_val})")
    
    # Calculate the final result to show the components
    final_calc_comp1 = sympy.S(1)/2 * p_gh_val
    final_calc_comp2 = sympy.S(1)/2 * p_ck_val
    print(f"P_MainEntrance = {final_calc_comp1} + {final_calc_comp2}")

    # Extract numerator and denominator for the final answer
    num, den = p_me_val.as_numer_denom()

    print(f"\nThus, the exact probability of reaching the Treasure Room from the Main Entrance is {p_me_val}, or as an irreducible fraction: {num}/{den}")

if __name__ == "__main__":
    solve_bran_castle_challenge()