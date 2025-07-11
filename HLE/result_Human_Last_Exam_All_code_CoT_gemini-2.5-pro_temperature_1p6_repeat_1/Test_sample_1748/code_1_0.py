def solve_rope_riddle():
    """
    Solves the magical rope riddle by iterating through possibilities.
    """
    initial_length = 60
    target_length = 15
    shrink_factor = 0.75

    # Iterate through possible lengths for the longer piece from the first cut.
    # It must be longer than half the total length.
    for p1a in range(initial_length // 2 + 1, initial_length):
        p1b = initial_length - p1a

        # After the cut, the longer piece (p1a) shrinks.
        # For the result to be an integer, p1a must be a multiple of 4.
        if p1a % 4 != 0:
            continue
        
        p1a_shrunk = p1a * shrink_factor

        # The piece used for the second cut (L2) must be > 30cm (i.e., > 2 * target_length)
        # to ensure the 15cm piece is the shorter one.
        # This piece must be p1a_shrunk, as p1b will always be < 30cm.
        if p1a_shrunk > 2 * target_length:
            rope_for_cut2 = p1a_shrunk
            p2b = rope_for_cut2 - target_length # The other piece from the second cut

            # For all final lengths to be integers, p2b (the longer part of the
            # second cut) must shrink to an integer. So, it must be a multiple of 4.
            if p2b > target_length and p2b % 4 == 0:
                # We found the unique solution. Now, we print the steps.
                p2b_shrunk = p2b * shrink_factor
                
                print("Here is the step-by-step solution:")
                print("-" * 40)
                print(f"1. Start with the initial rope of {initial_length}cm.")
                print(f"2. The first cut must be made at {p1a}cm, creating one piece of {p1a}cm and another of {p1b}cm.")
                print(f"3. The longer piece ({p1a}cm) shrinks by 25%.")
                print(f"   Equation: {p1a} * {shrink_factor} = {int(p1a_shrunk)}")
                print(f"4. The two pieces after the first cut are now {int(p1a_shrunk)}cm and {p1b}cm.")
                print("-" * 40)
                print("To verify, let's complete the second cut:")
                print(f"5. Take the {int(p1a_shrunk)}cm piece and cut it to get the {target_length}cm target piece.")
                print(f"6. The other piece will be {int(p1a_shrunk)} - {target_length} = {int(p2b)}cm.")
                print(f"7. This second piece ({int(p2b)}cm) is longer, so it shrinks by 25%.")
                print(f"   Equation: {int(p2b)} * {shrink_factor} = {int(p2b_shrunk)}")
                print("8. The process successfully yields a 15cm piece, and all final lengths (8cm, 15cm, 18cm) are integers.")
                print("-" * 40)
                
                answer = int(p1a_shrunk)
                print("The question asks for the length of the longer portion after the first cut.")
                print(f"This is the length of the shrunken piece from the first cut, which is {answer}cm.")

                return answer

solve_rope_riddle()