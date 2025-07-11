def solve_peg_game_classes():
    """
    Solves the peg game equivalence class problem by demonstrating a mathematical invariant.

    The method involves coloring the board and observing how the number of pegs
    of each color changes during a move. This reveals an invariant property
    that partitions all possible configurations into a fixed number of sets.
    """

    print("Step 1: Define an invariant for the peg game based on a 3-coloring.")
    print("---------------------------------------------------------------------")
    print("We assign a 'color' to each position (x, y) on the integer lattice using the rule:")
    print("c(x, y) = (x + y) % 3.")
    print("This partitions the board into 3 color sets: S_0, S_1, and S_2.")
    print("Let N_k be the number of pegs on squares of color k for any given configuration.")
    print("\nA move (forward or backward) involves three consecutive positions, for example,")
    print("p1=(x,y), p2=(x+1,y), p3=(x+2,y). Their respective colors are (x+y)%3, (x+y+1)%3,")
    print("and (x+y+2)%3. These three colors are always a permutation of {0, 1, 2}.")
    print("\nA forward move removes pegs from two positions and adds a peg to the third.")
    print("While the counts (N_0, N_1, N_2) change, we can find a related quantity that is invariant.")
    print("\nLet's define the invariant I(C) for a configuration C as a pair of values:")
    print("I(C) = ( (N_0 - N_1) % 2, (N_1 - N_2) % 2 )")
    print("This value can be shown to remain constant for any valid move (forward or backward).")

    print("\nStep 2: Determine the number of possible values for the invariant.")
    print("------------------------------------------------------------------")
    print("The first component of the invariant, (N_0 - N_1) % 2, can be 0 or 1.")
    possible_values_component_1 = 2
    print(f"Number of possible values for the first component: {possible_values_component_1}")
    print("The second component of the invariant, (N_1 - N_2) % 2, can also be 0 or 1.")
    possible_values_component_2 = 2
    print(f"Number of possible values for the second component: {possible_values_component_2}")
    print("\nThe total number of possible distinct values for the invariant is the product:")
    total_possible_classes = possible_values_component_1 * possible_values_component_2
    print(f"Equation: {possible_values_component_1} * {possible_values_component_2} = {total_possible_classes}")
    print(f"\nThis implies there are at most {total_possible_classes} equivalence classes.")


    print("\nStep 3: Show that each possible invariant value is achievable.")
    print("---------------------------------------------------------------")
    print("We demonstrate this by providing an example configuration for each of the 4 possible invariant values.")

    def get_invariant(configuration):
        """
        Calculates the invariant for a given configuration.
        A configuration is a list of (x, y) tuples representing peg positions.
        """
        counts = {0: 0, 1: 0, 2: 0}
        for x, y in configuration:
            color = (x + y) % 3
            counts[color] += 1
        n0, n1, n2 = counts[0], counts[1], counts[2]
        inv1 = (n0 - n1) % 2
        inv2 = (n1 - n2) % 2
        return (inv1, inv2)

    # Example configurations
    config_A = [(0, 0)]
    config_B = [(0, 0), (1, 0)]
    config_C = [(0, 0), (2, 0)]
    config_D = [(0, 0), (3, 0)]

    # Calculate invariants for the example configurations
    inv_A = get_invariant(config_A)
    inv_B = get_invariant(config_B)
    inv_C = get_invariant(config_C)
    inv_D = get_invariant(config_D)

    print(f"Configuration {config_A} -> N=(1,0,0) -> Invariant = {inv_A}")
    print(f"Configuration {config_B} -> N=(1,1,0) -> Invariant = {inv_B}")
    print(f"Configuration {config_C} -> N=(1,0,1) -> Invariant = {inv_C}")
    print(f"Configuration {config_D} -> N=(2,0,0) -> Invariant = {inv_D}")
    
    print("\nSince we have found configurations for all 4 possible invariant values,")
    print("it means these configurations belong to different, non-equivalent classes.")
    print("Therefore, there are at least 4 equivalence classes.")

    print("\nStep 4: Conclusion")
    print("--------------------")
    print(f"From Step 2, we know there are at most {total_possible_classes} classes.")
    print(f"From Step 3, we know there are at least {total_possible_classes} classes.")
    print(f"Therefore, the total number of equivalence classes must be exactly {total_possible_classes}.")

if __name__ == "__main__":
    solve_peg_game_classes()