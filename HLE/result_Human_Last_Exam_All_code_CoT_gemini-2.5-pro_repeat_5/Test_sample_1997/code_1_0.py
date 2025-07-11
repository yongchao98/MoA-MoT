def solve_genji_ko_chapter_39():
    """
    This function determines and prints the partition for Chapter 39 of The Tale of Genji.
    The incense pattern (Genji-kō) for Chapter 39, "Yūgiri", connects the incenses
    1 and 3, and separately connects 4 and 5, leaving 2 by itself.
    """
    
    # Define the partition as a list of lists.
    # Each inner list represents a subset of the partition.
    partition = [[1, 3], [2], [4, 5]]
    
    # Sort the subsets based on their smallest element to ensure increasing order.
    # The list is already sorted this way, but this step makes the logic explicit.
    partition.sort(key=lambda subset: min(subset))
    
    # Format each subset into a string like "{1, 3}"
    # The numbers within each subset are also sorted.
    subset_strings = []
    for subset in partition:
        # Build the string for a single subset, e.g., "{" + "1, 3" + "}"
        subset_str = "{" + ", ".join(map(str, sorted(subset))) + "}"
        subset_strings.append(subset_str)
        
    # Join the formatted subset strings and wrap them in outer braces
    # to represent the set of sets.
    final_output_string = "{" + ", ".join(subset_strings) + "}"
    
    print(final_output_string)

solve_genji_ko_chapter_39()