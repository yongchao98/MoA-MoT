def solve_rope_puzzle():
    """
    This function solves the magical rope puzzle by following a logical deduction process.
    """
    initial_rope_length = 60
    shrink_factor = 0.75
    target_length = 15

    print("Searching for the solution by analyzing the cuts...")
    print("-" * 30)

    # We iterate through possible lengths of the longer piece from the first cut.
    # It must be longer than half, so > 30.
    # It must be a multiple of 4 for the shrunk length to be an integer.
    for p1_long in range(32, initial_rope_length, 4):
        p1_short = initial_rope_length - p1_long
        
        # This is the piece that would be taken for the second cut.
        # It's the longer piece from the first cut, after shrinking.
        L1 = p1_long * shrink_factor

        # All lengths must be integers.
        if L1 != int(L1):
            continue
        L1 = int(L1)

        # From our analysis, the piece for the second cut (L1) must be > 30.
        if L1 <= 30:
            continue

        # And (L1 - 15) must be the longer piece of the second cut, and a multiple of 4.
        p2_long = L1 - target_length
        if p2_long > target_length and p2_long % 4 == 0:
            p2_shrunk = int(p2_long * shrink_factor)
            
            print("Solution Found!")
            print(f"1. First Cut: The {initial_rope_length}cm rope is cut into {p1_long}cm and {p1_short}cm.")
            
            # The question asks for the length of the longer portion AFTER the cut.
            answer = L1
            print(f"2. Magic Shrink: The longer piece ({p1_long}cm) shrinks by 25%.")
            
            # Outputting the numbers in the final equation as requested
            print(f"   Equation: {p1_long} * (1 - 25/100) = {answer}")
            
            print(f"3. Pieces after First Cut: We now have a {answer}cm rope and a {p1_short}cm rope.")
            print(f"4. Second Cut: The {answer}cm rope is cut into {p2_long}cm and {target_length}cm.")
            print(f"5. Magic Shrink: The longer piece ({p2_long}cm) shrinks to {p2_shrunk}cm.")
            print(f"6. Final Result: We successfully obtained a {target_length}cm piece of rope.")
            print("-" * 30)
            print(f"The length of the longer portion after the first cut is: {answer}cm")
            return

# Run the solver
solve_rope_puzzle()