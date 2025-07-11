def find_next_triplets():
    """
    Finds the next three triplets in the sequence based on a discovered mathematical rule.
    """
    # The rule is (a*x + b*y + c*z) % m == d
    # The letters X, Y, Z are represented by their 0-indexed position in the alphabet.
    a, b, c, m, d = 4, 10, 1, 14, 7

    # The last given triplet starts with 'N', so the next ones start with 'O'.
    x_char = 'O'
    x = ord(x_char) - ord('A')

    # Find all valid triplets starting with 'O'
    found_triplets = []
    for y_val in range(26):
        for z_val in range(26):
            # Check if the triplet (x, y_val, z_val) satisfies the rule
            if (a * x + b * y_val + c * z_val) % m == d:
                y_char = chr(ord('A') + y_val)
                z_char = chr(ord('A') + z_val)
                found_triplets.append(f"{x_char}{y_char}{z_char}")

    # The sequence is sorted alphabetically, so we sort the results
    found_triplets.sort()

    # Get the next three triplets
    next_three = found_triplets[:3]
    
    # Print the equation coefficients as requested
    print(f"The rule is based on the equation: ({a}*x + {b}*y + {c}*z) % {m} == {d}")
    
    # Print the final result
    print("The next three capital letters in the sequence are:")
    for triplet in next_three:
        print(triplet)

find_next_triplets()