def solve_genji_partition():
    """
    This function determines and prints the incense pattern partition for Chapter 39
    of the Tale of Genji.

    In the Genji-kō pattern for Chapter 39 ("Yūgiri"), the vertical bars
    for sticks 3 and 5 are connected, while bars 1, 2, and 4 are independent.
    This corresponds to a partition of the set {1, 2, 3, 4, 5}.
    """
    
    # Define the partition as a list of lists.
    # Each inner list is a subset of the partition.
    # {1}, {2}, {4}, {3, 5}
    partition_lists = [[1], [2], [4], [3, 5]]

    # Sort the numbers within each subset to ensure consistent ordering.
    for subset in partition_lists:
        subset.sort()

    # Sort the list of subsets based on their smallest element.
    # This fulfills the "sorted increasingly" requirement.
    sorted_partition = sorted(partition_lists, key=lambda s: s[0])

    # Format the sorted partition into the specified string format: {{a}, {b, c}, ...}
    subset_strings = []
    for subset in sorted_partition:
        # Convert each number in the subset to a string and join with ", "
        formatted_subset = "{" + ", ".join(map(str, subset)) + "}"
        subset_strings.append(formatted_subset)
    
    # Join all the formatted subsets with ", " and wrap with outer braces.
    final_output = "{" + ", ".join(subset_strings) + "}"
    
    print(final_output)

solve_genji_partition()