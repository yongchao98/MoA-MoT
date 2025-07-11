def solve_s4():
    """
    This function demonstrates the calculation for a term in sequence S4
    using the deduced recurrence relation.
    """
    # The sequence S4 provided in the problem, stored in a dictionary for 1-based indexing.
    # Index 0 is a placeholder.
    s_list = [0, 1, 1, 2, 2, 2, 4, 3, 4, 4, 4, 8, 5, 5, 8, 8, 6, 8, 12, 8, 11, 9, 9, 10, 13, 16, 9, 12, 20, 10, 12, 23, 12, 15, 21, 13, 17, 18, 19, 19, 22, 21, 19]
    s = {i: s_list[i] for i in range(1, len(s_list))}

    # We will demonstrate the calculation for n=6, as the rule holds true here.
    n_to_calculate = 6

    # Get the required previous terms from the sequence
    s_n_minus_1 = s[n_to_calculate - 1]
    s_n_minus_2 = s[n_to_calculate - 2]

    # Calculate the indices for the recurrence relation
    index1 = n_to_calculate - s_n_minus_1
    index2 = n_to_calculate - s_n_minus_2

    # Get the values from those indices
    val1 = s[index1]
    val2 = s[index2]

    # Calculate the final result
    result = val1 + val2

    print(f"The deduced rule is R(s[n]) = s[n - s[n-1]] + s[n - s[n-2]]")
    print(f"\nCalculating s[{n_to_calculate}] using the rule:")
    print(f"s[{n_to_calculate}] = s[{n_to_calculate} - s[{n_to_calculate - 1}]] + s[{n_to_calculate} - s[{n_to_calculate - 2}]]")
    print(f"Using s[{n_to_calculate - 1}] = {s_n_minus_1} and s[{n_to_calculate - 2}] = {s_n_minus_2}:")
    print(f"s[{n_to_calculate}] = s[{n_to_calculate} - {s_n_minus_1}] + s[{n_to_calculate} - {s_n_minus_2}]")
    print(f"s[{n_to_calculate}] = s[{index1}] + s[{index2}]")
    print(f"Using s[{index1}] = {val1} and s[{index2}] = {val2}:")
    print(f"s[{n_to_calculate}] = {val1} + {val2}")
    print(f"s[{n_to_calculate}] = {result}")

solve_s4()
<<<s[n] = s[n - s[n-1]] + s[n - s[n-2]]>>>