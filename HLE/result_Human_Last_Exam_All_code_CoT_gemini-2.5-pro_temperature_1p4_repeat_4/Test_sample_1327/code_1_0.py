def solve_sequence():
    """
    This function solves the number sequence puzzle by identifying a two-level pattern.
    """
    x = [2, 11, 23, 51, 119]

    # Step 1: Find the secondary sequence 'k' where x_n = 2 * x_{n-1} + k_{n-2} for n >= 3
    k = []
    # k_1 corresponds to the calculation of x_3
    k.append(x[2] - 2 * x[1])
    # k_2 corresponds to the calculation of x_4
    k.append(x[3] - 2 * x[2])
    # k_3 corresponds to the calculation of x_5
    k.append(x[4] - 2 * x[3])

    # print(f"The primary sequence is: {x}")
    # print(f"The secondary sequence 'k' is: {k}")

    # Step 2: Find the pattern in 'k'. We found k_m = 3 * k_{m-1} + 2
    # Calculate the next term in 'k'
    next_k = 3 * k[-1] + 2

    # Step 3: Calculate the next term in the primary sequence 'x'
    next_x = 2 * x[-1] + next_k

    # Step 4: Output the final equation with all numbers
    print(f"{x[0]}, {x[1]}, {x[2]}, {x[3]}, {x[4]}, ({next_x})")
    print("\nDerivation:")
    print(f"The pattern is x_n = 2 * x_(n-1) + k_(n-2) for n>=3, where k is a helper sequence.")
    print(f"{x[2]} = 2 * {x[1]} + {k[0]}")
    print(f"{x[3]} = 2 * {x[2]} + {k[1]}")
    print(f"{x[4]} = 2 * {x[3]} + {k[2]}")
    print("\nThe pattern for the helper sequence 'k' is k_m = 3 * k_(m-1) + 2.")
    print(f"{k[1]} = 3 * {k[0]} + 2")
    print(f"{k[2]} = 3 * {k[1]} + 2")
    print(f"So the next k is 3 * {k[2]} + 2 = {next_k}")
    print(f"Therefore, the next number in the main sequence is 2 * {x[4]} + {next_k} = {next_x}")

solve_sequence()