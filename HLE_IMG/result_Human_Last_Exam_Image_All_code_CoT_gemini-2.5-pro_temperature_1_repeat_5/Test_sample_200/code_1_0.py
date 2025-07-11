def solve_geology_problem():
    """
    Analyzes the geological events that could cause the observed rock transition.

    The transition from Rock A (likely deep-water shale) to Rock B (likely shallower-water sediment)
    represents a regression, which is caused by a relative fall in sea level.

    Let's evaluate the given options:
    i. Onset of glaciation: Locks up ocean water in ice sheets, causing global sea level to fall. This leads to regression. (Plausible)
    ii. Coastal subsidence: The land sinks, causing a relative sea level rise. This leads to transgression. (Not plausible)
    iii. Global warming: Melts ice and causes thermal expansion of water, leading to a global sea level rise and transgression. (Not plausible)
    iv. Decreased Mid-Ocean Ridge Activity: Ridges cool and contract, increasing ocean basin volume, causing global sea level to fall. This leads to regression. (Plausible)
    """

    plausible_events = {
        'i': "Onset of glaciation causes sea-level fall (regression).",
        'ii': "Coastal subsidence causes relative sea-level rise (transgression).",
        'iii': "Global warming causes sea-level rise (transgression).",
        'iv': "Decreased Mid-Ocean Ridge Activity causes sea-level fall (regression)."
    }

    correct_options = []
    for option, explanation in plausible_events.items():
        if "regression" in explanation:
            correct_options.append(option)
    
    # Sort the options for consistent output, though not strictly necessary
    correct_options.sort()

    # Format the answer as a comma-separated string with no spaces
    answer = ",".join(correct_options)
    print(f"The observed transition is a regression (shallowing-upward sequence).")
    print(f"The events causing a regression are:")
    for option in correct_options:
        print(f"- {option}: {plausible_events[option]}")
    
    print("\nFinal Answer:")
    print(answer)

solve_geology_problem()
<<<i,iv>>>