def find_next_number():
    """
    This function solves the number sequence puzzle by finding the pattern in the differences.
    """
    # The given sequence
    sequence = [2, 11, 23, 51, 119]
    print(f"The original sequence is: {', '.join(map(str, sequence))}")

    # Step 1: Calculate the differences between consecutive terms
    differences = [sequence[i+1] - sequence[i] for i in range(len(sequence)-1)]
    print(f"The differences between consecutive terms are: {', '.join(map(str, differences))}")
    print("Let's find the pattern in this new sequence of differences: 9, 12, 28, 68.")

    # Step 2: Analyze the pattern in the differences
    # d(n) = 2 * d(n-1) + c
    # 28 = 2 * 12 + 4
    # 68 = 2 * 28 + 12
    c1 = differences[2] - 2 * differences[1]
    c2 = differences[3] - 2 * differences[2]
    print(f"The pattern for the differences seems to be d(n) = 2 * d(n-1) + c(n).")
    print(f"Let's check:")
    print(f"{differences[2]} = 2 * {differences[1]} + {c1}")
    print(f"{differences[3]} = 2 * {differences[2]} + {c2}")
    
    # The added terms c(n) form a geometric progression: 4, 12, ...
    # The common ratio is 12 / 4 = 3
    common_ratio = c2 / c1
    print(f"The added terms are {c1}, {c2}, ... This is a geometric progression with a common ratio of {int(common_ratio)}.")

    # Step 3: Calculate the next added term
    next_c = c2 * common_ratio
    print(f"The next added term will be {c2} * {int(common_ratio)} = {int(next_c)}.")

    # Step 4: Calculate the next difference
    last_difference = differences[-1]
    next_difference = 2 * last_difference + next_c
    print(f"The next difference in the sequence will be 2 * {last_difference} + {int(next_c)} = {int(next_difference)}.")

    # Step 5: Calculate the next number in the original sequence
    last_number = sequence[-1]
    next_number = last_number + next_difference
    print(f"So, the number that should be filled in the parentheses is the last number plus the new difference.")
    print(f"The final equation is: {last_number} + {int(next_difference)} = {int(next_number)}")

find_next_number()