def solve_sequence():
    """
    This function solves the number sequence puzzle by identifying the underlying pattern.
    The pattern is found to be x_{n+1} = 2 * x_n + C_n, where the sequence of addends C_n
    follows its own geometric progression.
    """
    # The given sequence
    s = [2, 11, 23, 51, 119]

    # 1. We assume the pattern is x_{n+1} = 2 * x_n + C_n and find the addends C.
    # The pattern appears to start consistently from the second term.
    c1 = s[2] - 2 * s[1]  # 23 - 2*11 = 1
    c2 = s[3] - 2 * s[2]  # 51 - 2*23 = 5
    c3 = s[4] - 2 * s[3]  # 119 - 2*51 = 17
    addends = [c1, c2, c3]

    # 2. We find the pattern in the 'addends' sequence [1, 5, 17].
    # We look at the differences.
    diff1 = addends[1] - addends[0]  # 5 - 1 = 4
    diff2 = addends[2] - addends[1]  # 17 - 5 = 12

    # 3. The differences [4, 12] show a geometric progression with a ratio of 3.
    ratio = diff2 / diff1

    # 4. We predict the next difference and the next addend.
    next_diff = diff2 * ratio
    next_addend = addends[-1] + next_diff

    # 5. We calculate the next number in the original sequence.
    last_known_number = s[-1]
    next_number = 2 * last_known_number + next_addend

    # 6. Print the final equation as requested.
    print(f"The next term in the sequence is calculated as follows:")
    print(f"2 * {int(last_known_number)} + {int(next_addend)} = {int(next_number)}")

solve_sequence()