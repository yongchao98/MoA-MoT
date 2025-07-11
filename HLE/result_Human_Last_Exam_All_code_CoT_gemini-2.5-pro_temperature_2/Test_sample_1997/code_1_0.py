def solve_genji_ko():
    """
    Determines and prints the partition for Chapter 39 of the Tale of Genji.
    """
    # The partition is derived from the incense pattern for Chapter 39 ("Yugiri").
    # In this pattern, sticks 1 and 2 are linked, 3 and 4 are linked, and 5 is alone.
    partition = [{1, 2}, {3, 4}, {5}]

    # Sort the numbers within each set and convert them to lists for consistent ordering
    sorted_sets_as_lists = [sorted(list(s)) for s in partition]

    # Sort the list of lists based on the first element of each inner list
    sorted_partition = sorted(sorted_sets_as_lists)

    # Format the output string as a set of sets
    # e.g., {{1,2},{3,4},{5}}
    output_parts = []
    for part in sorted_partition:
        # Create string for each inner set, e.g., "{1,2}"
        inner_set_str = "{" + ",".join(map(str, part)) + "}"
        output_parts.append(inner_set_str)
    
    final_output = "{" + ",".join(output_parts) + "}"
    
    print(final_output)

solve_genji_ko()