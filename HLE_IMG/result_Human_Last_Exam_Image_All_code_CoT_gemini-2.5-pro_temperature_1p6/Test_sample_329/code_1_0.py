def solve_simulation_puzzle():
    """
    This function analyzes the four triple-slit simulations to find the unique parameter
    for each, and then combines their codes into a single four-digit sequence.

    The parameters are:
    w: slit width (inversely related to beam spread)
    u: wave number/phase (affects interference and tilt)
    h: height (related to intensity/brightness at the source)

    The encoding is:
    w1: 1, w2: 2, w3: 3
    u1: 4, u2: 5, u3: 6
    h1: 7, h2: 8, h3: 9
    """

    # Analysis and results for each simulation
    # Simulation 1:
    # - Identical: Slit intensities are equal (h1=h2=h3).
    # - Monotonic: Beam spread decreases left-to-right, so width w increases (w1<w2<w3).
    # - Unique: By elimination, 'u' is unique. The central interference is disrupted,
    #           which implies the central slit's phase is different (u1=u3 != u2).
    # Unique parameter is u2.
    code_1 = 5

    # Simulation 2:
    # - Identical: Beam spreads are equal (w1=w2=w3).
    # - Monotonic: Slit intensities decrease left-to-right (h1>h2>h3).
    # - Unique: By elimination, 'u' is unique. The pattern is asymmetric,
    #           favoring the left side, implying the left slit's phase is different (u2=u3 != u1).
    # Unique parameter is u1.
    code_2 = 4

    # Simulation 3:
    # - Identical: Slit intensities are equal (h1=h2=h3).
    # - Monotonic: The entire pattern is strongly tilted, implying a monotonic
    #              change in phase (u1<u2<u3 or vice versa).
    # - Unique: By elimination, 'w' is unique. The central beam looks different from the
    #           outer two, implying the central slit's width is different (w1=w3 != w2).
    # Unique parameter is w2.
    code_3 = 2

    # Simulation 4:
    # - Identical: The pattern is symmetric and not tilted, so phase is constant (u1=u2=u3).
    # - Monotonic: Beam spread decreases left-to-right, so width w increases (w1<w2<w3).
    # - Unique: By elimination, 'h' is unique. The leftmost slit is dimmer
    #           than the other two (h2=h3 != h1).
    # Unique parameter is h1.
    code_4 = 7

    # Combine the codes into a single sequence
    final_sequence = f"{code_1}{code_2}{code_3}{code_4}"

    print(f"Code for Simulation 1 (Unique parameter u2): {code_1}")
    print(f"Code for Simulation 2 (Unique parameter u1): {code_2}")
    print(f"Code for Simulation 3 (Unique parameter w2): {code_3}")
    print(f"Code for Simulation 4 (Unique parameter h1): {code_4}")
    print(f"\nThe combined four-digit sequence is: {final_sequence}")

solve_simulation_puzzle()
<<<5427>>>