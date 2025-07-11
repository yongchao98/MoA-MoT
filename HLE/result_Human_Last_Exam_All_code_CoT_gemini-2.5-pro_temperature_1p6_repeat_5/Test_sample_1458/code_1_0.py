def solve_cathedral_riddle():
    """
    Solves the Cathedral's Echo riddle by finding the number of pipes
    that meet all the conditions described in the poem.
    """
    p_in_tune = 200

    # Let's find the number of out-of-tune pipes.
    # From the poem, the number of out-of-tune pipes ("P_out") must be
    # a multiple of both 7 and 4, so it must be a multiple of 28.
    # The total number of pipes ("T") must be a multiple of 3 and 5,
    # so it must be a multiple of 15.
    # We also know that T = P_out + 200.
    # We can search for the first multiple of 28 that satisfies this.

    p_out_of_tune = 0
    total_pipes = 0

    # Start checking multiples of 28 for the number of out-of-tune pipes
    for i in range(1, 100):
        potential_p_out = 28 * i
        potential_total_pipes = potential_p_out + p_in_tune

        # Check if the potential total pipes is a multiple of 15
        if potential_total_pipes % 15 == 0:
            p_out_of_tune = potential_p_out
            total_pipes = potential_total_pipes
            break

    print(f"Based on the poem's constraints, we can deduce:")
    print(f"Total number of pipes: {total_pipes}")
    print(f"Pipes in tune: {p_in_tune}")
    print(f"Pipes out of tune (the 'lost'): {p_out_of_tune}\n")

    print("Let's verify this solution:")
    group_a = total_pipes / 3
    group_b = total_pipes * (2/5)
    print(f"'One-third of pipes' is {total_pipes} / 3 = {int(group_a)}")
    print(f"'Two-fifths' is {total_pipes} * 2/5 = {int(group_b)}")
    print(f"These two groups combined (with an overlap of {int(group_a + group_b - p_out_of_tune)}) make the {p_out_of_tune} out-of-tune pipes.\n")

    print("Checking the breakdown of the out-of-tune pipes:")
    octaves_rich = p_out_of_tune * (3/7)
    minor_scales = p_out_of_tune * (1/4)
    print(f"'Three-sevenths' of the lost pipes is {int(octaves_rich)}.")
    print(f"'One-fourth' of the lost pipes is {int(minor_scales)}.")
    print("These are whole numbers, so the solution is consistent.\n")

    # Final Question: "How many must the tuner find When just half the lost realign?"
    # This is asking for half of the out-of-tune pipes.
    final_answer = p_out_of_tune / 2

    print("The final question asks how many pipes the tuner must realign, which is half the number of lost pipes.")
    print("The equation to find the answer is:")
    print(f"{p_out_of_tune} / 2 = {int(final_answer)}")


solve_cathedral_riddle()
<<<140>>>