def solve_mancala_puzzle():
    """
    Solves the Mancala score difference problem by using a parity argument.
    """
    # Initial game state
    p1_pits = [0, 2, 0, 0, 2, 0]
    p1_store = 22
    p2_pits = [1, 0, 0, 0, 0, 0]
    p2_store = 21

    # Step 1: Calculate the total number of stones.
    total_p1_pits = sum(p1_pits)
    total_p2_pits = sum(p2_pits)
    total_stones = total_p1_pits + p1_store + total_p2_pits + p2_store

    print("Step 1: Calculate the total number of stones in the game.")
    print(f"Player 1's pits: {total_p1_pits} stones")
    print(f"Player 1's store: {p1_store} stones")
    print(f"Player 2's pits: {total_p2_pits} stones")
    print(f"Player 2's store: {p2_store} stones")
    print("---")
    print("The total number of stones is the sum of all stones on the board.")
    print(f"Total Stones = {total_p1_pits} + {p1_store} + {total_p2_pits} + {p2_store} = {total_stones}")
    print("-" * 30)

    # Step 2: Apply the parity argument.
    print("Step 2: Apply the parity argument.")
    print("Let the final scores for Player 1 and Player 2 be S1 and S2.")
    print("At the end of the game, all stones are in the stores, so S1 + S2 must equal the total number of stones.")
    print(f"Therefore, S1 + S2 = {total_stones}.")

    is_total_even = (total_stones % 2 == 0)
    print(f"\nThe total number of stones, {total_stones}, is an {'even' if is_total_even else 'odd'} number.")
    print("For the sum of two integers (S1 and S2) to be even, they must have the same parity (i.e., both are even, or both are odd).")
    print("The difference between two integers of the same parity is always an even number.")
    print("Therefore, the final score difference, |S1 - S2|, must be an even number.")
    print("-" * 30)

    # Step 3: Evaluate the answer choices.
    print("Step 3: Evaluate the answer choices.")
    choices = {"A": 0, "B": 1, "C": 2, "D": 3, "E": 4, "F": 5}
    unobtainable_choices_count = 0

    for choice, value in choices.items():
        if value % 2 != 0:
            print(f"- Score difference {value} (Choice {choice}) is an odd number, so it is unobtainable.")
            unobtainable_choices_count += 1
        else:
            print(f"- Score difference {value} (Choice {choice}) is an even number, so it is potentially obtainable.")
    print("-" * 30)

    # Step 4: Formulate the conclusion.
    print("Step 4: Formulate the conclusion.")
    if unobtainable_choices_count > 1:
        print("Since the score differences of 1, 3, and 5 are all unobtainable, more than one of the listed score differences is unobtainable.")
    else:
        # This case won't be reached based on the logic, but is included for completeness.
        print("There is only one unobtainable score difference in the list.")

solve_mancala_puzzle()