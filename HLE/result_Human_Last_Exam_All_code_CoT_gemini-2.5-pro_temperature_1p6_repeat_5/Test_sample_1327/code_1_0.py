def solve_sequence():
    """
    This function analyzes the sequence to find the next term.
    """
    # The given sequence
    s = [2, 11, 23, 51, 119]

    print("Analyzing the sequence:", s)
    print("Let's test the hypothesis that each term is related to the previous term by a formula like: Term(n) = 3 * Term(n-1) - C")
    print("-" * 30)

    # We establish the pattern starting from the 3rd term
    # Term 3: 23
    c2 = 3 * s[1] - s[2]
    print(f"For the 3rd term ({s[2]}): {s[2]} = 3 * {s[1]} - {c2}")

    # Term 4: 51
    c3 = 3 * s[2] - s[3]
    print(f"For the 4th term ({s[3]}): {s[3]} = 3 * {s[2]} - {c3}")

    # Term 5: 119
    c4 = 3 * s[3] - s[4]
    print(f"For the 5th term ({s[4]}): {s[4]} = 3 * {s[3]} - {c4}")

    # Now we have a sequence of the subtracted numbers C
    c_sequence = [c2, c3, c4]
    print("-" * 30)
    print("This reveals a sequence of subtracted numbers:", c_sequence)

    # Let's find the pattern in this new sequence by looking at the differences
    diff1 = c_sequence[1] - c_sequence[0]
    diff2 = c_sequence[2] - c_sequence[1]
    print(f"The differences in this subtraction sequence are {diff1} and {diff2}.")
    
    # The differences are 8 and 16. It seems the difference doubles each time.
    print("The difference doubles each time (8 * 2 = 16).")
    next_diff = diff2 * 2
    print(f"The next difference should be {diff2} * 2 = {next_diff}.")

    # Now we can find the next number in the subtraction sequence.
    next_c = c_sequence[2] + next_diff
    print(f"So, the next number in the subtraction sequence is {c_sequence[2]} + {next_diff} = {next_c}.")
    print("-" * 30)

    # Finally, we can calculate the next number in the original sequence.
    last_known_term = s[4]
    next_term = 3 * last_known_term - next_c
    print("Therefore, the next number in the original sequence is found by applying the same rule:")
    print(f"Next Term = 3 * {last_known_term} - {next_c} = {3 * last_known_term} - {next_c} = {next_term}")

solve_sequence()