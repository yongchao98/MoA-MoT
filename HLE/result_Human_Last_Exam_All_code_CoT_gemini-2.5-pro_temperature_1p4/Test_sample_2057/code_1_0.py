def solve_hat_puzzle():
    """
    This function prints the solution to the hat puzzle.
    """

    # The number of people who deduce their hat color in Round 3.
    # This corresponds to the number of people with white hats.
    num_yes = 4

    # One possible distribution of hats that fits the conditions.
    # Other valid solutions are rotations of this sequence.
    # B = Black, W = White
    hat_distribution = "BWBWBWBBW"

    print(f"How many people replied 'Yes.'?")
    print(num_yes)
    print("\nHow are black and white hats distributed around the table?")
    print("One possible configuration is:")
    
    # We represent the circle of people and hats
    people = [f"Person {i+1}" for i in range(9)]
    hats = list(hat_distribution)
    
    # Create a visual representation of the table
    output_lines = []
    for i in range(9):
        hat_color = "Black" if hats[i] == 'B' else "White"
        output_lines.append(f"{people[i]}: {hat_color}")
        
    print("\n".join(output_lines))
    
    # As per the logic, the 4 people with White hats are the ones who say "Yes".
    # All 4 are in a symmetric position (a BWB triplet), so they solve it simultaneously.
    
solve_hat_puzzle()