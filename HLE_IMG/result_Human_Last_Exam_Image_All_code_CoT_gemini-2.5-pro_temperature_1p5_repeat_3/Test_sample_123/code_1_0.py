import sympy

def solve_bran_castle_escape():
    """
    Calculates the exact probability of reaching the Treasure Room from the Main Entrance
    using absorbing Markov chain analysis.
    """
    # 1. Define the states
    # Transient states are the rooms where you can still move around.
    # Absorbing states are the final outcomes (Treasure or Doom).
    transient_states = ['ME', 'GH', 'CK', 'SP', 'Lib', 'KH', 'TC', 'QR']
    absorbing_states = ['TR', 'VL']

    # 2. Set up the transition matrices Q and R using exact fractions.
    # Q: Transitions between transient states.
    # R: Transitions from transient to absorbing states.
    # Rows/Cols for Q: ME, GH, CK, SP, Lib, KH, TC, QR
    # Cols for R: TR, VL
    
    # Using sympy.Rational for exact fraction arithmetic
    Q = sympy.Matrix([
        [0, sympy.Rational(1, 2), sympy.Rational(1, 2), 0, 0, 0, 0, 0],                            # ME -> GH(1/2), CK(1/2)
        [sympy.Rational(1, 4), 0, 0, sympy.Rational(1, 4), sympy.Rational(1, 4), sympy.Rational(1, 4), 0, 0], # GH -> ME(1/4), SP(1/4), Lib(1/4), KH(1/4)
        [sympy.Rational(1, 2), 0, 0, 0, 0, sympy.Rational(1, 2), 0, 0],                            # CK -> ME(1/2), KH(1/2)
        [0, sympy.Rational(1, 4), 0, 0, 0, 0, sympy.Rational(1, 4), 0],                            # SP -> GH(1/4), TC(1/4)
        [0, sympy.Rational(1, 3), 0, 0, 0, sympy.Rational(1, 3), 0, sympy.Rational(1, 3)],          # Lib -> GH(1/3), KH(1/3), QR(1/3)
        [0, sympy.Rational(1, 4), sympy.Rational(1, 4), 0, sympy.Rational(1, 4), 0, 0, sympy.Rational(1, 4)], # KH -> GH(1/4), CK(1/4), Lib(1/4), QR(1/4)
        [0, 0, 0, sympy.Rational(1, 2), 0, 0, 0, 0],                                               # TC -> SP(1/2)
        [0, 0, 0, 0, sympy.Rational(1, 3), sympy.Rational(1, 3), 0, 0]                             # QR -> Lib(1/3), KH(1/3)
    ])

    R = sympy.Matrix([
        [0, 0],                               # ME -> None
        [0, 0],                               # GH -> None
        [0, 0],                               # CK -> None
        [sympy.Rational(1, 4), sympy.Rational(1, 4)],  # SP -> TR(1/4), VL(1/4)
        [0, 0],                               # Lib -> None
        [0, 0],                               # KH -> None
        [0, sympy.Rational(1, 2)],                     # TC -> VL(1/2)
        [sympy.Rational(1, 3), 0]                      # QR -> TR(1/3)
    ])

    # 3. Perform the matrix calculations.
    # Get the number of transient states
    num_transient = Q.shape[0]

    # Create the identity matrix I
    I = sympy.eye(num_transient)

    # Calculate the fundamental matrix N = (I - Q)^-1
    N = (I - Q).inv()

    # Calculate the absorption probability matrix B = N * R
    B = N * R

    # 4. Extract and print the result.
    # The starting state is Main Entrance (ME), which is index 0.
    start_state_index = transient_states.index('ME')
    
    # The target success state is Treasure Room (TR), which is index 0 in the absorbing list.
    treasure_room_index = absorbing_states.index('TR')
    
    # Get the probability from the B matrix
    probability_fraction = B[start_state_index, treasure_room_index]
    
    numerator = probability_fraction.p
    denominator = probability_fraction.q

    print("This problem is solved by calculating the absorption probability in a Markov chain.")
    print("The system of rooms forms a state graph, and we want to find the probability of ending in the 'Treasure Room' state.")
    print("\nUsing the matrix method for absorbing Markov chains, we compute B = (I - Q)^-1 * R.")
    print("The resulting probability for starting at the Main Entrance is:")
    
    # Output the final equation with the numbers filled in
    print("\nP(Treasure Room | Start = Main Entrance) =", numerator, "/", denominator)


solve_bran_castle_escape()