import sympy

def solve_bran_castle_challenge():
    """
    Solves the Bran Castle escape challenge by setting up and solving a system of linear equations.
    """
    # Define the 8 unknown probabilities for the non-absorbing states
    # x0: Main Entrance, x1: Great Hall, x2: Castle Kitchen, x3: Secret Passage,
    # x4: Library, x5: Knights' Hall, x6: Torture Chamber, x7: Queen's Room
    x = sympy.symbols('x0:8')
    room_names = [
        "P(Main Entrance)", "P(Great Hall)", "P(Castle Kitchen)", "P(Secret Passage)",
        "P(Library)", "P(Knights' Hall)", "P(Torture Chamber)", "P(Queen's Room)"
    ]

    # System of equations based on the principle:
    # P(current) = sum(prob_to_next * P(next))
    # P(Treasure Room) = 1
    # P(Vampire's Lair) = 0
    # Equations are written in the form: expression = 0

    # Main Entrance (2 exits: Great Hall, Castle Kitchen)
    # x0 = 1/2*x1 + 1/2*x2
    eq0 = x[0] - sympy.S(1)/2 * x[1] - sympy.S(1)/2 * x[2]

    # Great Hall (4 exits: Main Entrance, Secret Passage, Library, Knights' Hall)
    # x1 = 1/4*x0 + 1/4*x3 + 1/4*x4 + 1/4*x5
    eq1 = x[1] - sympy.S(1)/4 * x[0] - sympy.S(1)/4 * x[3] - sympy.S(1)/4 * x[4] - sympy.S(1)/4 * x[5]

    # Castle Kitchen (2 exits: Main Entrance, Knights' Hall)
    # x2 = 1/2*x0 + 1/2*x5
    eq2 = x[2] - sympy.S(1)/2 * x[0] - sympy.S(1)/2 * x[5]

    # Secret Passage (4 exits: Great Hall, Torture Chamber, Treasure Room, Vampire's Lair)
    # x3 = 1/4*x1 + 1/4*x6 + 1/4*1 + 1/4*0
    eq3 = x[3] - sympy.S(1)/4 * x[1] - sympy.S(1)/4 * x[6] - sympy.S(1)/4

    # Library (3 exits: Great Hall, Queen's Room, Knights' Hall)
    # x4 = 1/3*x1 + 1/3*x7 + 1/3*x5
    eq4 = x[4] - sympy.S(1)/3 * x[1] - sympy.S(1)/3 * x[7] - sympy.S(1)/3 * x[5]

    # Knights' Hall (4 exits: Great Hall, Library, Queen's Room, Castle Kitchen)
    # x5 = 1/4*x1 + 1/4*x4 + 1/4*x7 + 1/4*x2
    eq5 = x[5] - sympy.S(1)/4 * x[1] - sympy.S(1)/4 * x[4] - sympy.S(1)/4 * x[7] - sympy.S(1)/4 * x[2]

    # Torture Chamber (2 exits: Secret Passage, Vampire's Lair)
    # x6 = 1/2*x3 + 1/2*0
    eq6 = x[6] - sympy.S(1)/2 * x[3]

    # Queen's Room (3 exits: Library, Knights' Hall, Treasure Room)
    # x7 = 1/3*x4 + 1/3*x5 + 1/3*1
    eq7 = x[7] - sympy.S(1)/3 * x[4] - sympy.S(1)/3 * x[5] - sympy.S(1)/3

    print("Step 1: Define the system of linear equations for the probability of success from each room.")
    print("Let P(Room) be the probability of reaching the Treasure Room from that room.\n")
    print("The equations are:")
    print(f"{room_names[0]} = 1/2 * {room_names[1]} + 1/2 * {room_names[2]}")
    print(f"{room_names[1]} = 1/4 * {room_names[0]} + 1/4 * {room_names[3]} + 1/4 * {room_names[4]} + 1/4 * {room_names[5]}")
    print(f"{room_names[2]} = 1/2 * {room_names[0]} + 1/2 * {room_names[5]}")
    print(f"{room_names[3]} = 1/4 * {room_names[1]} + 1/4 * {room_names[6]} + 1/4 * 1 + 1/4 * 0")
    print(f"{room_names[4]} = 1/3 * {room_names[1]} + 1/3 * {room_names[7]} + 1/3 * {room_names[5]}")
    print(f"{room_names[5]} = 1/4 * {room_names[1]} + 1/4 * {room_names[4]} + 1/4 * {room_names[7]} + 1/4 * {room_names[2]}")
    print(f"{room_names[6]} = 1/2 * {room_names[3]} + 1/2 * 0")
    print(f"{room_names[7]} = 1/3 * {room_names[4]} + 1/3 * {room_names[5]} + 1/3 * 1")
    
    print("\nStep 2: Solve the system of equations to find the probabilities.")
    equations = [eq0, eq1, eq2, eq3, eq4, eq5, eq6, eq7]
    solution = sympy.solve(equations, x)
    
    prob_me = solution[x[0]]
    numerator = prob_me.p
    denominator = prob_me.q
    
    print("\nStep 3: The solution for the probability of starting from the Main Entrance is:")
    print(f"{room_names[0]} = {numerator} / {denominator}")

solve_bran_castle_challenge()
<<<109/161>>>