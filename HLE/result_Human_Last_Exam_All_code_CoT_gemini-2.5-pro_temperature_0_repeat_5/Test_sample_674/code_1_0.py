def solve_sequence():
    """
    This function generates the sequence S4 based on its recurrence relation
    and prints the result.
    """
    s = {1: 1, 2: 1}
    n_terms = 42

    # Generate the sequence up to n_terms
    for n in range(3, n_terms + 1):
        try:
            # The deduced rule for the sequence S4
            # s[n] = s[s[n-2]] + s[n-s[n-1]]
            s[n] = s[s[n - 2]] + s[n - s[n - 1]]
        except KeyError as e:
            print(f"Error: Index {e} not found in sequence while calculating s[{n}].")
            return

    # Print the full sequence
    sequence_str = ", ".join(str(s[i]) for i in range(1, n_terms + 1))
    print(f"S4 = {sequence_str}")

    # Print the equation for the last term as an example
    n = n_terms
    s_n_minus_1 = s[n - 1]
    s_n_minus_2 = s[n - 2]
    index1 = s[n - 2]
    index2 = n - s[n - 1]
    val1 = s[index1]
    val2 = s[index2]
    result = s[n]

    print("\nRule: s[n] = s[s[n-2]] + s[n-s[n-1]]")
    print(f"For n = {n}:")
    print(f"s[{n}] = s[s[{n-2}]] + s[{n}-s[{n-1}]]")
    print(f"s[{n}] = s[{s_n_minus_2}] + s[{n}-{s_n_minus_1}]")
    print(f"s[{n}] = s[{index1}] + s[{index2}]")
    print(f"s[{n}] = {val1} + {val2}")
    print(f"s[{n}] = {result}")

solve_sequence()
<<<R(s[n]) = s[s[n-2]] + s[n-s[n-1]]>>>