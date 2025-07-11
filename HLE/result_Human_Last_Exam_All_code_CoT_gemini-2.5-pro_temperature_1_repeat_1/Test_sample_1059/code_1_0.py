def solve_closepact_problem():
    """
    Solves the topological problem by identifying which sets are necessarily 'closepact'.

    The logic is as follows:
    1. A set is 'closepact in itself' if it's an H-closed space.
    2. For the given examples (subspaces of metric spaces), H-closed is equivalent to compact.
    3. For subsets of R or C, compact is equivalent to closed and bounded.
    4. Each option is evaluated based on whether it is necessarily compact.
    """

    # Dictionary to store the evaluation of each choice.
    # True if the set is necessarily compact, False otherwise.
    choices = {
        'A': False,  # R is not bounded.
        'B': False,  # Z is not bounded.
        'C': True,   # Any finite set is compact.
        'D': False,  # Not closed (limit point 0 is not in the set).
        'E': False,  # Not necessarily compact (e.g., sequence converging to a point not in the set).
        'F': False,  # Not necessarily closed (if limit point is not included).
        'G': True,   # Bounded, and adding the limit point makes it closed. Hence, compact.
        'H': True,   # A convergent sequence is bounded, and adding the limit makes it closed. Hence, compact.
        'I': False,  # An open interval is not closed.
        'J': True,   # A closed interval is closed and bounded, the definition of compact in R.
        'K': False,  # Not necessarily closed (e.g., an open interval).
        'L': False,  # A non-measurable set in R is not closed, hence not compact.
        'M': True    # The Cantor set is closed and bounded, hence compact.
    }

    # Construct the final answer string by concatenating the letters of the correct choices.
    # The letters are processed in alphabetical order to ensure a consistent output.
    result_string = ""
    for letter in sorted(choices.keys()):
        if choices[letter]:
            result_string += letter

    # Print the final result string.
    print(result_string)

solve_closepact_problem()