def solve_sequence():
    """
    This function solves the number sequence puzzle.
    The sequence is: 2, 11, 23, 51, 119, ...

    Step-by-step analysis:
    1. The pattern is assumed to be a_{n+1} = 2 * a_n + c_n.
    2. We calculate the sequence c_n starting from the second term.
    3. We find the pattern within the sequence c_n to predict its next term.
    4. Finally, we calculate the next term of the original sequence.
    """
    a = [2, 11, 23, 51, 119]

    # Step 1: Find the sequence c_n for n >= 2
    # The relationship is a_{n+1} = 2 * a_n + c_{n+1}
    # We will analyze the pattern from the third term of 'a'
    # a[2] = 2*a[1] + c_2 => 23 = 2*11 + c_2
    c = []
    for i in range(1, len(a) - 1):
        c_val = a[i+1] - 2 * a[i]
        c.append(c_val)

    # c is now [1, 5, 17]

    # Step 2: Find the pattern in c
    # The differences are:
    diffs_c = []
    for i in range(len(c) - 1):
        diff = c[i+1] - c[i]
        diffs_c.append(diff)
    
    # diffs_c is [4, 12]

    # Step 3: Assume the differences form a geometric progression
    ratio = diffs_c[1] / diffs_c[0]  # 12 / 4 = 3
    next_diff = diffs_c[-1] * ratio

    # Step 4: Calculate the next term in c
    next_c = c[-1] + next_diff

    # Step 5: Calculate the next term in a
    next_a = 2 * a[-1] + next_c
    
    print("The pattern is that each term is twice the previous term plus an increasing number.")
    print("Sequence: 2, 11, 23, 51, 119, ...")
    print("23 = 2 * 11 + 1")
    print("51 = 2 * 23 + 5")
    print("119 = 2 * 51 + 17")
    print("The added numbers are 1, 5, 17. The differences between these are 4 and 12.")
    print("The differences (4, 12) form a geometric progression with a ratio of 3.")
    print(f"The next difference is 12 * 3 = {12*3}.")
    print(f"The next number to add is 17 + 36 = {17+36}.")
    print(f"Therefore, the next term in the sequence is 2 * 119 + 53.")
    print()
    print("The final equation is:")
    # Using integer casting to ensure the output is clean
    print(f"{2} * {a[-1]} + {int(next_c)} = {int(next_a)}")

solve_sequence()
<<<A>>>