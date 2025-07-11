def get_genji_ko_partition():
    """
    This function determines and formats the Genji-kō partition for Chapter 39.

    The incense pattern for Chapter 39, "Yūgiri", has three groups:
    1. The first stick is alone: {1}
    2. The second stick is alone: {2}
    3. The third, fourth, and fifth sticks are connected: {3, 4, 5}

    The problem specifies a particular ordering for the sets in the final output,
    based on the example for chapter 20. We will follow that ordering rule.
    """
    # The raw partition for Chapter 39.
    partition = [{1}, {2}, {3, 4, 5}]

    # Create a mapping from each number to the set it belongs to.
    num_to_set = {num: s for s in partition for num in s}

    ordered_partition = []
    seen_numbers = set()

    # Iterate from 1 to 5 to establish the correct order of sets
    # as per the problem's implicit rule.
    for i in range(1, 6):
        if i not in seen_numbers:
            # Find the set containing the current number i.
            current_set = num_to_set[i]
            # Add this set to our ordered result.
            ordered_partition.append(current_set)
            # Mark all numbers in this set as seen to avoid re-adding.
            seen_numbers.update(current_set)

    # Format the ordered partition into the final string representation.
    # e.g., "{{1}, {2}, {3, 4, 5}}"
    # Within each set, numbers are sorted for consistent output.
    formatted_sets = []
    for s in ordered_partition:
        inner_str = "{" + ", ".join(map(str, sorted(list(s)))) + "}"
        formatted_sets.append(inner_str)

    final_output_string = "{" + ", ".join(formatted_sets) + "}"

    print(final_output_string)

get_genji_ko_partition()