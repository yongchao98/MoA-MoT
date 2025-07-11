def solve_sequence():
    """
    This function solves the number sequence puzzle by identifying a two-level recursive pattern.
    """
    # The given sequence
    x = [2, 11, 23, 51, 119]

    # Step 1: Assume the pattern is x_{n+1} = 2 * x_n + c_n and find the constants 'c'.
    # We observe the rule seems to apply consistently from the second term onwards.
    # c2 is the constant needed to get from x2=11 to x3=23
    # c3 is the constant needed to get from x3=23 to x4=51
    # c4 is the constant needed to get from x4=51 to x5=119
    c2 = x[2] - 2 * x[1]  # 23 - 2 * 11 = 1
    c3 = x[3] - 2 * x[2]  # 51 - 2 * 23 = 5
    c4 = x[4] - 2 * x[3]  # 119 - 2 * 51 = 17

    constants = [c2, c3, c4]

    # Step 2: Find the pattern in the sequence of constants: 1, 5, 17.
    # We hypothesize a recursive pattern c_{i+1} = 3 * c_i + 2.
    # Let's verify:
    # 3 * c2 + 2 = 3 * 1 + 2 = 5 (which is c3)
    # 3 * c3 + 2 = 3 * 5 + 2 = 17 (which is c4)
    # The pattern holds true.

    # Step 3: Calculate the next constant, c5.
    c5 = 3 * c4 + 2

    # Step 4: Calculate the next term in the main sequence, x6.
    x6 = 2 * x[4] + c5
    
    # Step 5: Print the results including the final equation.
    print(f"The original sequence is: {', '.join(map(str, x))}, ...")
    print(f"The pattern discovered is x(n+1) = 2 * x(n) + c(n) for n>=2.")
    print(f"The sequence of constants 'c' is: {', '.join(map(str, constants))}, ...")
    print(f"The pattern for constants is c(n+1) = 3 * c(n) + 2.")
    print(f"The next constant is c5 = 3 * {c4} + 2 = {c5}.")
    print(f"Therefore, the next term in the sequence is calculated as:")
    print(f"2 * {x[4]} + {c5} = {x6}")

solve_sequence()