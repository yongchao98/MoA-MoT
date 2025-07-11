def solve_sequence():
    """
    This function solves the number sequence puzzle.
    The pattern is a_n = 3 * a_{n-1} - c_{n-1}, where c follows its own pattern.
    """
    # The given sequence
    sequence = [2, 11, 23, 51, 119]
    print(f"Original sequence: {', '.join(map(str, sequence))}, ()")
    print("-" * 30)

    # Let's find the pattern.
    # We observe that a_n = 3 * a_{n-1} - c
    # Let's find the sequence of 'c'
    c_values = []
    print("Step 1: Discover the underlying pattern.")
    # The pattern starts from the 3rd term in the sequence
    # For a_3 = 23
    c1 = 3 * sequence[1] - sequence[2]
    c_values.append(c1)
    print(f"{sequence[2]} = 3 * {sequence[1]} - {c1}")

    # For a_4 = 51
    c2 = 3 * sequence[2] - sequence[3]
    c_values.append(c2)
    print(f"{sequence[3]} = 3 * {sequence[2]} - {c2}")

    # For a_5 = 119
    c3 = 3 * sequence[3] - sequence[4]
    c_values.append(c3)
    print(f"{sequence[4]} = 3 * {sequence[3]} - {c3}")
    print(f"The sequence of subtracted numbers is: {c_values}")
    print("-" * 30)

    # Now let's find the pattern in the c_values sequence
    print("Step 2: Analyze the sequence of subtracted numbers.")
    diff1 = c_values[1] - c_values[0]
    print(f"The difference between the first two subtracted numbers is: {c_values[1]} - {c_values[0]} = {diff1}")
    diff2 = c_values[2] - c_values[1]
    print(f"The difference between the next two is: {c_values[2]} - {c_values[1]} = {diff2}")
    print("The difference doubles each time (8, 16, ...).")
    print("-" * 30)
    
    # Calculate the next difference and the next c value
    print("Step 3: Predict the next terms.")
    next_diff = diff2 * 2
    print(f"The next difference will be: {diff2} * 2 = {next_diff}")
    next_c = c_values[-1] + next_diff
    print(f"The next subtracted number will be: {c_values[-1]} + {next_diff} = {next_c}")
    print("-" * 30)

    # Calculate the next number in the original sequence
    print("Step 4: Calculate the final missing number in the sequence.")
    next_term = 3 * sequence[-1] - next_c
    print(f"Next term = 3 * {sequence[-1]} - {next_c} = {3 * sequence[-1]} - {next_c} = {next_term}")
    
    # Print the final completed sequence
    print("-" * 30)
    final_sequence = sequence + [next_term]
    print(f"The final sequence is: {', '.join(map(str, sequence))}, ({next_term})")

solve_sequence()
<<<A>>>