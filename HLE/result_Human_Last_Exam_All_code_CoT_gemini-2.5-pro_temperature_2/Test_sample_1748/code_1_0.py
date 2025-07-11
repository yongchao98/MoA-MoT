def solve_rope_puzzle():
    """
    This function solves the rope puzzle by iterating through possible first cuts
    and checking if they satisfy all the problem's constraints.
    """
    initial_length = 60
    # The first cut produces a shorter piece 'a' and a longer piece 'b'.
    # We iterate through possible integer lengths for the longer piece 'b'.
    # 'b' must be longer than half the total length.
    for b in range(initial_length // 2 + 1, initial_length):
        
        # Constraint: After the longer piece 'b' shrinks, its length must be an integer.
        # b * 0.75 = b * 3/4 must be an integer, which means b must be a multiple of 4.
        if b % 4 != 0:
            continue

        a = initial_length - b
        shrunken_b = int(b * 0.75) # This is the length of the longer portion after the first cut.

        # Let's call the rope for the second cut L1.
        # Based on logic, we must use the shrunken longer piece for the second cut.
        L1 = shrunken_b
        
        # Constraint for the second cut:
        # To get a 15cm piece, it must be the shorter part.
        # So, the other part (L1 - 15) must be longer than 15. This means L1 > 30.
        if L1 <= 30:
            continue

        # Constraint: The longer piece from the second cut, (L1 - 15), must also
        # result in an integer length after shrinking.
        # (L1 - 15) * 0.75 must be an integer, so (L1 - 15) must be a multiple of 4.
        if (L1 - 15) % 4 == 0:
            
            # We found the solution. Let's print the steps.
            print("The solution is found through the following steps:")
            print(f"1. The initial 60cm rope is cut into two pieces: {a}cm and {b}cm.")
            print(f"2. The longer portion ({b}cm) shrinks by 25%.")
            print(f"Equation: {b} * (1 - 0.25) = {shrunken_b}")
            print(f"3. The length of the longer portion after this first cut is {shrunken_b}cm.")
            
            # We continue the steps to verify this is a valid solution.
            print(f"4. This {shrunken_b}cm piece is then cut to get the 15cm rope.")
            piece_y = L1 - 15
            print(f"5. The cut produces a {15}cm piece and a {piece_y}cm piece, successfully obtaining the goal.")
            print("\nFinal Answer:")
            print(f"The length of the longer portion after the first cut is {shrunken_b}cm.")
            
            # The puzzle is solved, so we can exit the loop.
            return shrunken_b

# Run the solver and get the final answer.
final_answer = solve_rope_puzzle()
# The value is returned, but the print inside the function handles the output.
# No need to print final_answer here as the function does it.