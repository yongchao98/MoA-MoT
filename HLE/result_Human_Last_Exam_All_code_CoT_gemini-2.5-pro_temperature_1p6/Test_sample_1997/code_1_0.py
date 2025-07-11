def solve_genji_partition():
    """
    This function determines and prints the Genji-mon partition for Chapter 39.

    The Genji-mon for Chapter 39 ("Yūgiri") connects the middle three sticks,
    leaving the first and last separate. This corresponds to the partition
    of the set {1, 2, 3, 4, 5} into the subsets {1}, {2, 3, 4}, and {5}.
    """
    
    # Define the partition for Chapter 39 as a list of sets.
    # Using frozenset to make them hashable if needed, though list works fine here.
    partition = [
        {1},
        {2, 3, 4},
        {5}
    ]

    # Sort the list of sets based on the minimum element in each set.
    # This fulfills the "sorted increasingly" requirement.
    sorted_partition = sorted(partition, key=min)

    # Prepare the string representation for each subset.
    subset_strings = []
    for subset in sorted_partition:
        # Sort elements within each subset for consistent formatting.
        sorted_elements = sorted(list(subset))
        # Format as "{e1,e2,...}"
        subset_strings.append("{" + ", ".join(map(str, sorted_elements)) + "}")

    # Combine the subset strings into the final output format "{ {s1}, {s2}, ... }"
    final_output_string = "{" + ", ".join(subset_strings) + "}"

    print(final_output_string)

solve_genji_partition()