def solve_genji_partition():
    """
    This function determines and prints the partition for Chapter 39 of "The Tale of Genji".

    In the Genji-kō incense patterns, the 5 sticks are numbered 1 to 5 from right to left.
    For Chapter 39 ("Yūgiri"), the pattern connects sticks 2 and 3, while sticks 1, 4,
    and 5 are independent. This corresponds to the partition: {1}, {2, 3}, {4}, {5}.

    The code formats this partition as a set of sets, sorted increasingly by the
    minimum element of each subset.
    """

    # The partition for Chapter 39, with subsets already sorted by their minimum element.
    partition_ch39 = [[1], [2, 3], [4], [5]]

    # This list will hold the string representation of each subset, e.g., "{1}", "{2, 3}"
    subset_strings = []

    for subset in partition_ch39:
        # Sort the numbers within each subset for consistent ordering, e.g., [3, 2] -> [2, 3]
        subset.sort()
        # Convert the list of numbers to a list of strings
        string_numbers = [str(num) for num in subset]
        # Format the subset as a string, e.g., "{2, 3}"
        subset_str = "{" + ", ".join(string_numbers) + "}"
        subset_strings.append(subset_str)

    # Combine all subset strings into the final representation, e.g., "{{1}, {2, 3}, ...}"
    final_output = "{" + ", ".join(subset_strings) + "}"

    print(final_output)

solve_genji_partition()