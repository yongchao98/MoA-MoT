def find_closepact_sets():
    """
    Identifies which of the given sets are necessarily "closepact".

    The problem defines "closepact" as what is commonly known as H-closed.
    For the given choices, which are all subspaces of metric spaces, this
    property is equivalent to being compact. For subsets of R or C, a set is
    compact if and only if it is closed and bounded (Heine-Borel Theorem).

    This function encodes the analysis for each choice.
    """

    # A dictionary to hold the reasoning for each choice.
    # True means the set is necessarily compact, False otherwise.
    analysis = {
        'A': False,  # The set of real numbers is not bounded.
        'B': False,  # The set of integers is not bounded.
        'C': True,   # A finite set in a metric space is always closed and bounded.
        'D': False,  # Not closed, as it's missing the limit point 0.
        'E': False,  # Not necessarily compact (e.g., a sequence in Q converging to an irrational number).
        'F': False,  # Not necessarily closed, as the limit point may be missing.
        'G': True,   # A bounded monotonic sequence converges. The set including the limit is closed and bounded.
        'H': True,   # A convergent sequence is bounded. The set including the limit is closed and bounded.
        'I': False,  # An open interval is not closed.
        'J': True,   # A closed interval is closed and bounded by definition.
        'K': False,  # Not necessarily closed (e.g., (0,1) is bounded and measurable but not closed).
        'L': False,  # Not closed (a compact set must be closed, and all closed sets are measurable).
        'M': True    # The Cantor set is famously closed and bounded.
    }

    # Collect the letters of the choices that are necessarily compact.
    correct_choices = []
    for choice in sorted(analysis.keys()):
        if analysis[choice]:
            correct_choices.append(choice)

    # Print the final answer as a single string.
    print("".join(correct_choices))

find_closepact_sets()