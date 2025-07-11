def solve_sequence():
    """
    This function solves the number sequence puzzle by finding a recursive pattern.
    """
    # The given sequence
    s = [2, 11, 23, 51, 119]
    print(f"The sequence is: {', '.join(map(str, s))}, ...")

    # The pattern is of the form: next_num = 2 * prev_num + addend
    # Let's find the addends. We'll start from the second term.
    # The calculation for the first term (11) is often an exception or a seed.
    addends = []
    # Starting from s[2] = 23 = 2 * s[1] + c1
    for i in range(1, len(s) - 1):
        addend = s[i+1] - 2 * s[i]
        addends.append(addend)
        print(f"{s[i+1]} = 2 * {s[i]} + {addend}")

    # Now we analyze the pattern in the addends: [1, 5, 17]
    print(f"\nThe sequence of addends is: {', '.join(map(str, addends))}")

    # Let's find the pattern in the differences of the addends
    addend_diffs = []
    for i in range(len(addends) - 1):
        diff = addends[i+1] - addends[i]
        addend_diffs.append(diff)
    print(f"The differences of the addends are: {', '.join(map(str, addend_diffs))}")

    # The differences [4, 12] form a geometric progression with ratio 3
    ratio = addend_diffs[1] / addend_diffs[0]
    print(f"This is a geometric progression with a common ratio of {int(ratio)}.")

    # Calculate the next difference
    next_diff = addend_diffs[-1] * ratio
    print(f"The next difference will be {addend_diffs[-1]} * {int(ratio)} = {int(next_diff)}.")

    # Calculate the next addend
    next_addend = addends[-1] + next_diff
    print(f"The next addend will be {addends[-1]} + {int(next_diff)} = {int(next_addend)}.")

    # Calculate the final number in the original sequence
    last_num = s[-1]
    next_num = 2 * last_num + next_addend

    print("\nTherefore, the final equation to find the next number is:")
    # Here we output each number in the final equation as requested.
    print(f"2 * {last_num} + {int(next_addend)} = {int(next_num)}")

solve_sequence()
<<<A>>>