def solve_sequence():
    """
    Solves the number sequence puzzle by identifying a recursive pattern.
    """
    # The given sequence
    x = [2, 11, 23, 51, 119]
    print(f"The given sequence is: {', '.join(map(str, x))}")

    # Step 1: Propose the rule x_n = 3 * x_{n-1} - k
    # and find the values for k.
    k_1 = 3 * x[2-1] - x[3-1] # k for calculating x_3
    k_2 = 3 * x[3-1] - x[4-1] # k for calculating x_4
    k_3 = 3 * x[4-1] - x[5-1] # k for calculating x_5
    
    k_sequence = [k_1, k_2, k_3]
    print(f"The pattern can be described as x_n = 3 * x_(n-1) - k, where k changes.")
    print(f"Let's find the values of k:")
    print(f"For {x[2]}: {x[2]} = 3 * {x[1]} - {k_1}")
    print(f"For {x[3]}: {x[3]} = 3 * {x[2]} - {k_2}")
    print(f"For {x[4]}: {x[4]} = 3 * {x[3]} - {k_3}")
    print(f"The sequence of subtracted numbers 'k' is: {k_sequence}")

    # Step 2: Find the pattern in the sequence 'k'
    diff1 = k_sequence[1] - k_sequence[0]
    diff2 = k_sequence[2] - k_sequence[1]
    print(f"\nLet's analyze the sequence k = {k_sequence}.")
    print(f"The differences are {diff1} and {diff2}.")
    print("This is a geometric progression (each difference is 2 times the previous one).")

    # Step 3: Predict the next k
    next_diff = diff2 * 2
    next_k = k_sequence[-1] + next_diff
    print(f"The next difference should be {diff2} * 2 = {next_diff}.")
    print(f"So the next k is {k_sequence[-1]} + {next_diff} = {next_k}.")

    # Step 4: Calculate the next term in the original sequence
    next_x = 3 * x[-1] - next_k
    print(f"\nFinally, we calculate the next term in the original sequence:")
    print(f"Next Term = 3 * {x[-1]} - {next_k} = {3 * x[-1]} - {next_k} = {next_x}")
    print(f"The complete sequence is: {', '.join(map(str, x))}, ({next_x})")

solve_sequence()
<<<A>>>