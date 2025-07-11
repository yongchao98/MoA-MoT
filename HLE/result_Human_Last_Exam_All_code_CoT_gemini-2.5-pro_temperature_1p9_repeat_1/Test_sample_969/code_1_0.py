def solve_sequence():
    """
    Solves the number sequence puzzle based on the identified pattern.
    """
    # The given starting sequence (first 9 elements)
    given_sequence = [3, 2, 1, 2, 3, 3, 3, 2, 2]

    # The prefix that defines the pattern
    prefix = [3, 2, 1, 2]

    # Generate the full sequence based on the rule
    # The rule: Prefix, followed by blocks of each prefix element repeated 3 times.
    full_sequence = list(prefix)
    for k in prefix:
        full_sequence.extend([k] * 3)

    # Determine the next 4 elements
    start_index = len(given_sequence)
    next_four_elements = full_sequence[start_index : start_index + 4]

    # Combine the given sequence and the found elements
    final_sequence = given_sequence + next_four_elements

    # Print the full sequence as requested
    print("The complete sequence is:")
    # Using a list comprehension to convert numbers to strings for join
    print(" ".join([str(num) for num in final_sequence]))

solve_sequence()