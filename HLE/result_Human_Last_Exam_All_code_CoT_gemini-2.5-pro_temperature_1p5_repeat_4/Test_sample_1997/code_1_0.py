def solve_genji_ko_partition():
    """
    This function determines and prints the incense pattern for Chapter 39
    of The Tale of Genji as a partition of the set {1, 2, 3, 4, 5}.

    The Genji-kō (incense pattern) for Chapter 39, "Yūgiri", connects the first
    three sticks, leaving the fourth and fifth sticks independent.
    """

    # The partition is represented as a list of lists.
    # The outer list and inner lists are sorted to ensure the final
    # output is "sorted increasingly" as requested.
    # {1, 2, 3} are one group.
    # {4} is a separate group.
    # {5} is a separate group.
    partition_ch39 = [[1, 2, 3], [4], [5]]

    # Format each inner list into a string that looks like a set.
    # For example, [1, 2, 3] becomes "{1, 2, 3}".
    inner_set_strings = []
    for subset in partition_ch39:
        # The join function requires strings, so we map each integer to a string.
        subset_content = ", ".join(map(str, subset))
        inner_set_strings.append(f"{{{subset_content}}}")

    # Join the formatted inner set strings with a comma and space.
    # Then, wrap the entire result in curly braces to represent the outer set.
    final_output = f"{{{', '.join(inner_set_strings)}}}"

    print(final_output)

solve_genji_ko_partition()