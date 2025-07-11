def solve_rope_puzzle():
    """
    This function solves the rope puzzle by iterating through possible cuts
    and checking if they meet the problem's constraints.
    """
    initial_length = 60
    target_length = 15
    shrink_factor = 0.75

    # The first longer piece 'p1_longer' must be > 30.
    # For the length after the first cut to be an integer, 'p1_longer' must be a multiple of 4.
    # We start checking from the first multiple of 4 greater than 30, which is 32.
    for p1_longer in range(32, initial_length, 4):
        p1_shorter = initial_length - p1_longer

        # Calculate the total rope length after the first cut and shrink.
        length_after_cut1 = p1_shorter + p1_longer * shrink_factor

        # This length must be an integer, which is guaranteed by our loop step.
        length_after_cut1 = int(length_after_cut1)

        # For the second cut, the 15cm piece must be the shorter one.
        # This means the other piece must be longer.
        p2_longer = length_after_cut1 - target_length
        p2_shorter = target_length

        if p2_longer > p2_shorter:
            # Calculate the total rope length after the second cut and shrink.
            length_after_cut2 = p2_shorter + p2_longer * shrink_factor

            # Check if this final length is an integer.
            if length_after_cut2 == int(length_after_cut2):
                # We found the solution.
                # The question asks for the length of the longer portion after the first cut.
                longer_portion_after_cut1 = p1_longer * shrink_factor

                print("The problem can be solved with the following steps:")
                print(f"1. First cut: The 60cm rope is cut into a {p1_longer}cm piece and a {p1_shorter}cm piece.")
                print(f"2. After shrinking, the total rope length is {p1_shorter} + ({p1_longer} * {shrink_factor}) = {length_after_cut1}cm.")
                print(f"3. Second cut: The {length_after_cut1}cm rope is cut into a {p2_longer}cm piece and a {p2_shorter}cm piece.")
                print(f"4. We get the desired {p2_shorter}cm piece, and the final rope length is {int(length_after_cut2)}cm, which is an integer.")
                print("\nThe question is: What will be the length of the longer portion after the first cut?")
                print("This is the length of the first longer piece after it shrinks.")
                print("\nFinal Equation:")
                print(f"{p1_longer} * {shrink_factor} = {int(longer_portion_after_cut1)}")
                return

solve_rope_puzzle()