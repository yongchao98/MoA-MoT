def solve_genji_chapter_39():
    """
    This function determines the incense pattern (partition) for Chapter 39
    of The Tale of Genji and prints it in the specified format.
    """

    # The Genji-kō pattern for Chapter 39 ("Yūgiri") corresponds to the
    # partition where elements 3 and 5 are grouped, and 1, 2, and 4 are singletons.
    # We represent this partition as a list of sets.
    partition = [
        {1},
        {2},
        {4},
        {3, 5}
    ]

    # Sort the partition. The key for sorting is the minimum element in each subset.
    # This ensures an increasing order as specified in the prompt.
    sorted_partition = sorted(partition, key=min)

    # Format the output string to match the requested format, e.g., {{1},{2},{3,5},{4}}
    # First, format each inner set. We sort the elements within each set for a
    # consistent and clean representation (e.g., {5, 3} becomes {3,5}).
    inner_set_strings = []
    for subset in sorted_partition:
        # Sort numbers inside each subset before printing
        sorted_subset_list = sorted(list(subset))
        # Format as a string, e.g., [1] -> "{1}", [3, 5] -> "{3,5}"
        subset_str = '{' + ','.join(map(str, sorted_subset_list)) + '}'
        inner_set_strings.append(subset_str)

    # Join the formatted inner sets and wrap with outer braces
    final_output_string = '{' + ','.join(inner_set_strings) + '}'

    # Print the final result
    print(final_output_string)

solve_genji_chapter_39()