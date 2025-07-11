def solve_sequence():
    """
    This function solves the given number sequence puzzle.
    The sequence is 2, 11, 23, 51, 119, ( ).
    """
    # The original sequence
    sequence = [2, 11, 23, 51, 119]
    print(f"The original sequence is: {sequence}")

    # Step 1: Find the pattern. We assume a pattern of the form:
    # next_number = 2 * current_number + added_value
    # Let's find the sequence of 'added_value' starting from the second term.

    c1 = sequence[2] - 2 * sequence[1]  # 23 - 2 * 11 = 1
    c2 = sequence[3] - 2 * sequence[2]  # 51 - 2 * 23 = 5
    c3 = sequence[4] - 2 * sequence[3]  # 119 - 2 * 51 = 17
    added_sequence = [c1, c2, c3]
    
    print(f"\nThe relationship between terms seems to be `S_n = 2 * S_(n-1) + C_n`:")
    print(f"23 = 2 * 11 + {c1}")
    print(f"51 = 2 * 23 + {c2}")
    print(f"119 = 2 * 51 + {c3}")
    print(f"This gives a secondary sequence of added numbers: {added_sequence}")

    # Step 2: Find the pattern in the added_sequence: [1, 5, 17]
    d1 = added_sequence[1] - added_sequence[0]  # 5 - 1 = 4
    d2 = added_sequence[2] - added_sequence[1]  # 17 - 5 = 12
    differences = [d1, d2]
    
    print(f"\nThe differences in this secondary sequence are: {differences}")
    
    # The differences [4, 12] form a geometric progression with a ratio of 3.
    ratio = differences[1] / differences[0]
    print(f"This difference sequence is a geometric progression with a ratio of {int(ratio)}.")

    # Step 3: Predict the next terms based on the discovered pattern.
    next_difference = differences[1] * ratio
    print(f"The next difference will be {differences[1]} * {int(ratio)} = {int(next_difference)}.")
    
    next_added_value = added_sequence[2] + next_difference
    print(f"The next number in the added sequence will be {added_sequence[2]} + {int(next_difference)} = {int(next_added_value)}.")

    # Step 4: Calculate the final missing number in the original sequence.
    last_known_number = sequence[-1]
    final_answer = 2 * last_known_number + next_added_value
    
    print("\nFinally, we calculate the next number in the original sequence:")
    print(f"The final equation is: 2 * {last_known_number} + {int(next_added_value)} = {int(final_answer)}")

solve_sequence()