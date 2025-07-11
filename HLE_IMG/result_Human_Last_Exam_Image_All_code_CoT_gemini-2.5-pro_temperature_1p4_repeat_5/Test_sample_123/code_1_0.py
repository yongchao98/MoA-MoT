import sympy

def solve_bran_castle_challenge():
    """
    Solves the Bran Castle escape challenge by setting up and solving a system of
    linear equations representing the probabilities of reaching the treasure room
    from each location in the castle.
    """

    # Define the symbols for the probability of reaching the treasure from each room.
    # We use prefixes to make them easier to identify.
    # P0: Main Entrance, P1: Great Hall, P2: Castle Kitchen, P3: Secret Passage,
    # P4: Library, P5: Knights' Hall, P6: Queen's Room, P7: Torture Chamber
    P0, P1, P2, P3, P4, P5, P6, P7 = sympy.symbols('P0 P1 P2 P3 P4 P5 P6 P7')

    # The probabilities for the absorbing states are known constants.
    # P_TreasureRoom = 1
    # P_VampiresLair = 0

    # Set up the system of linear equations based on the connections and probabilities.
    # Each equation represents the probability of reaching the treasure from one room,
    # which is the sum of (probability to move to next room * probability from that next room).

    # Eq1: Main Entrance (connects to Great Hall, Castle Kitchen)
    # Exits: 2, Prob: 1/2 each
    eq1 = sympy.Eq(P0, sympy.Rational(1, 2) * P1 + sympy.Rational(1, 2) * P2)

    # Eq2: Great Hall (connects to Main Entrance, Secret Passage, Library, Knights' Hall)
    # Exits: 4, Prob: 1/4 each
    eq2 = sympy.Eq(P1, sympy.Rational(1, 4) * P0 + sympy.Rational(1, 4) * P3 + \
                       sympy.Rational(1, 4) * P4 + sympy.Rational(1, 4) * P5)

    # Eq3: Castle Kitchen (connects to Main Entrance, Knights' Hall)
    # Exits: 2, Prob: 1/2 each
    eq3 = sympy.Eq(P2, sympy.Rational(1, 2) * P0 + sympy.Rational(1, 2) * P5)

    # Eq4: Secret Passage (connects to Great Hall, Torture Chamber, Treasure Room, Vampire's Lair)
    # Exits: 4, Prob: 1/4 each. P_TreasureRoom=1, P_VampiresLair=0
    eq4 = sympy.Eq(P3, sympy.Rational(1, 4) * P1 + sympy.Rational(1, 4) * P7 + \
                       sympy.Rational(1, 4) * 1 + sympy.Rational(1, 4) * 0)

    # Eq5: Library (connects to Great Hall, Queen's Room, Knights' Hall)
    # Exits: 3, Prob: 1/3 each
    eq5 = sympy.Eq(P4, sympy.Rational(1, 3) * P1 + sympy.Rational(1, 3) * P6 + \
                       sympy.Rational(1, 3) * P5)

    # Eq6: Knights' Hall (connects to Great Hall, Library, Queen's Room, Castle Kitchen)
    # Exits: 4, Prob: 1/4 each
    eq6 = sympy.Eq(P5, sympy.Rational(1, 4) * P1 + sympy.Rational(1, 4) * P4 + \
                       sympy.Rational(1, 4) * P6 + sympy.Rational(1, 4) * P2)

    # Eq7: Queen's Room (connects to Library, Knights' Hall, Treasure Room)
    # Exits: 3, Prob: 1/3 each. P_TreasureRoom=1
    eq7 = sympy.Eq(P6, sympy.Rational(1, 3) * P4 + sympy.Rational(1, 3) * P5 + \
                       sympy.Rational(1, 3) * 1)

    # Eq8: Torture Chamber (connects to Secret Passage, Vampire's Lair)
    # Exits: 2, Prob: 1/2 each. P_VampiresLair=0
    eq8 = sympy.Eq(P7, sympy.Rational(1, 2) * P3 + sympy.Rational(1, 2) * 0)

    # Solve the system of equations for all unknown probabilities.
    solution = sympy.solve((eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8),
                             (P0, P1, P2, P3, P4, P5, P6, P7))

    # Extract the probabilities needed for the final calculation display.
    prob_main_entrance = solution[P0]
    prob_great_hall = solution[P1]
    prob_castle_kitchen = solution[P2]

    # Print the explanation and the final equation.
    print("The probability of reaching the Treasure Room from the Main Entrance is P_ME.")
    print("This probability is determined by the chances of success from the rooms it connects to:")
    print("the Great Hall (P_GH) and the Castle Kitchen (P_CK).")
    print("\nThe relationship is defined by the equation:")
    print("P_ME = 1/2 * P_GH + 1/2 * P_CK\n")

    print("By solving the complete system of equations for the castle, we find the following probabilities:")
    print(f"Probability from Great Hall (P_GH) = {prob_great_hall}")
    print(f"Probability from Castle Kitchen (P_CK) = {prob_castle_kitchen}\n")

    print("Substituting these values into the equation for the Main Entrance gives:")
    print(f"P_ME = 1/2 * ({prob_great_hall}) + 1/2 * ({prob_castle_kitchen})")
    
    # Show the intermediate calculation steps as fractions
    term1 = sympy.Rational(1, 2) * prob_great_hall
    term2 = sympy.Rational(1, 2) * prob_castle_kitchen
    print(f"P_ME = {term1} + {term2}")
    
    final_prob = prob_main_entrance
    print(f"P_ME = {final_prob}\n")

    print(f"The exact probability of reaching the Treasure Room from the Main Entrance is {final_prob}.")

if __name__ == '__main__':
    solve_bran_castle_challenge()