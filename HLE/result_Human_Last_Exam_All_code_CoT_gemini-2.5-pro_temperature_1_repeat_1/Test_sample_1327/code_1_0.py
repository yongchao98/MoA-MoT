def solve_sequence():
    """
    This function solves the number sequence puzzle by identifying a two-level pattern.
    """
    # The given sequence
    sequence = [2, 11, 23, 51, 119]
    print(f"The initial sequence is: {', '.join(map(str, sequence))}, ...")

    # The pattern is x_n = 2 * x_{n-1} + k, where k itself follows a pattern.
    # We establish this pattern starting from the third term.
    # For the term 23:
    k1 = sequence[2] - 2 * sequence[1]
    # For the term 51:
    k2 = sequence[3] - 2 * sequence[2]
    # For the term 119:
    k3 = sequence[4] - 2 * sequence[3]
    
    print("\nStep 1: Find the relationship between terms.")
    print(f"{sequence[2]} = 2 * {sequence[1]} + {k1}")
    print(f"{sequence[3]} = 2 * {sequence[2]} + {k2}")
    print(f"{sequence[4]} = 2 * {sequence[3]} + {k3}")
    
    k_sequence = [k1, k2, k3]
    print(f"\nStep 2: Find the pattern in the remainders {k_sequence}.")
    # The pattern for the k_sequence is k_new = 3 * k_old + 2
    print(f"We observe that {k_sequence[1]} = 3 * {k_sequence[0]} + 2")
    print(f"And {k_sequence[2]} = 3 * {k_sequence[1]} + 2")

    # Step 3: Calculate the next remainder.
    next_k = 3 * k_sequence[-1] + 2
    print(f"\nStep 3: Calculate the next remainder in the series.")
    print(f"The next remainder is 3 * {k_sequence[-1]} + 2 = {next_k}")

    # Step 4: Calculate the next term in the main sequence.
    last_term = sequence[-1]
    next_term_in_sequence = 2 * last_term + next_k
    print(f"\nStep 4: Calculate the final number in the sequence.")
    print(f"The final equation is: 2 * {last_term} + {next_k} = {next_term_in_sequence}")

solve_sequence()
<<<A>>>