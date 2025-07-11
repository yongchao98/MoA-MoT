def solve_sequence():
    """
    Calculates the S4 sequence based on the Hofstadter Q-sequence rule
    and prints the result.
    """
    s = {}
    try:
        # Base cases
        s[1] = 1
        s[2] = 1

        # Generate sequence up to n=42
        n = 42
        for i in range(3, n + 1):
            idx1 = i - s[i - 1]
            idx2 = i - s[i - 2]
            
            # Check if indices are valid before accessing
            if idx1 not in s or idx2 not in s:
                print(f"Error: Cannot compute s[{i}] due to invalid index.")
                return

            s[i] = s[idx1] + s[idx2]

        # Get values needed for the final equation
        s_n_minus_1 = s[n - 1]
        s_n_minus_2 = s[n - 2]
        term1_index = n - s_n_minus_1
        term2_index = n - s_n_minus_2
        term1_value = s[term1_index]
        term2_value = s[term2_index]
        final_value = s[n]

        # Print the sequence
        print("Generated S4 Sequence (Hofstadter Q-sequence):")
        sequence_str = ", ".join(str(s[i]) for i in range(1, n + 1))
        print(sequence_str)
        
        # Print the rule for the last term
        print("\nRule for the last term (s[42]):")
        print(f"s[{n}] = s[{n} - s[{n - 1}]] + s[{n} - s[{n - 2}]]")
        print(f"s[{n}] = s[{n} - {s_n_minus_1}] + s[{n} - {s_n_minus_2}]")
        print(f"s[{n}] = s[{term1_index}] + s[{term2_index}]")
        print(f"s[{n}] = {term1_value} + {term2_value} = {final_value}")

    except KeyError as e:
        print(f"Error computing sequence: A required key {e} was not found.")

solve_sequence()