def solve_geology_problem():
    """
    This function determines the plausible geological events for the observed rock transition.

    The transition from Rock A (deep-water deposit, likely shale) to Rock B
    (shallower-water deposit, likely sandstone/limestone) indicates a marine regression,
    which is a relative fall in sea level.

    We will evaluate which of the given options cause a sea-level fall.
    """

    # Dictionary of events and their effect on sea level.
    # 'fall' = regression, 'rise' = transgression.
    events = {
        'i': {'name': 'Onset of glaciation', 'effect': 'fall'},
        'ii': {'name': 'Coastal subsidence', 'effect': 'rise'},
        'iii': {'name': 'Global warming', 'effect': 'rise'},
        'iv': {'name': 'Decreased Mid-Ocean Ridge Activity', 'effect': 'fall'}
    }

    # The observed phenomenon is a regression (sea-level fall).
    observed_effect = 'fall'

    # Find all events that match the observed effect.
    plausible_events = []
    for key, value in events.items():
        if value['effect'] == observed_effect:
            plausible_events.append(key)

    # Sort the keys for consistent output.
    plausible_events.sort()

    # Format the answer as a comma-separated string with no spaces.
    answer = ",".join(plausible_events)
    print(answer)

solve_geology_problem()