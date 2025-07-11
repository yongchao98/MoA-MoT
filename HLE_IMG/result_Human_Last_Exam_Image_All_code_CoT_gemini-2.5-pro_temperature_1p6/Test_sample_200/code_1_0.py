def solve_geology_question():
    """
    Analyzes the geological events that could cause the observed rock transition.

    The transition from Rock A (deep-water sediment) to Rock B (shallower-water sediment)
    indicates a marine regression, or a fall in relative sea level.

    We evaluate the options:
    i. Onset of glaciation: Locks water in ice sheets, causing global sea level to fall. This causes a regression. Plausible.
    ii. Coastal subsidence: The land sinks, causing a relative sea level rise (transgression). Not plausible.
    iii. Global warming: Melts ice and expands ocean water, causing global sea level to rise (transgression). Not plausible.
    iv. Decreased Mid-Ocean Ridge Activity: Ridges cool and shrink, increasing ocean basin volume and causing global sea level to fall. This causes a regression. Plausible.

    The correct options are i and iv.
    """
    plausible_events = ['i', 'iv']
    answer = ",".join(plausible_events)
    print(answer)

solve_geology_question()