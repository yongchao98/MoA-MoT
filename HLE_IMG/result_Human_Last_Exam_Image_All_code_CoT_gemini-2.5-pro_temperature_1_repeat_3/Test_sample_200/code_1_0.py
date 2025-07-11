def solve_geology_problem():
    """
    This script analyzes the geological transition from Rock A to Rock B
    and determines the plausible causal events.
    """

    # Step 1: Interpret the geological evidence from the image.
    # The transition from Rock A (dark, fine-grained, deep-water deposit) to
    # Rock B (lighter, coarser, shallow-water deposit) indicates a
    # shallowing-upward sequence. This is called a marine regression, which is
    # caused by a relative fall in sea level.
    cause_of_sequence = "Relative sea-level fall (Regression)"

    # Step 2: Define the options and analyze their effect on sea level.
    options = {
        'i': {'event': 'Onset of glaciation', 'effect': 'Sea-level fall'},
        'ii': {'event': 'Coastal subsidence', 'effect': 'Relative sea-level rise'},
        'iii': {'event': 'Global warming', 'effect': 'Sea-level rise'},
        'iv': {'event': 'Decreased Mid-Ocean Ridge Activity', 'effect': 'Sea-level fall'}
    }

    # Step 3: Identify the options that match the required cause (sea-level fall).
    plausible_options = []
    for key, value in options.items():
        # A regression is caused by a sea-level fall.
        if "fall" in value['effect']:
            plausible_options.append(key)

    # Step 4: Format the final answer as a comma-separated string.
    # The list needs to be sorted for consistent output, although in this case it's already in order.
    plausible_options.sort()
    final_answer = ",".join(plausible_options)

    # Print the final answer
    print(final_answer)

solve_geology_problem()
<<<i,iv>>>