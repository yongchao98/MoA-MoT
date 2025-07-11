def solve_geology_problem():
    """
    Analyzes the geological events that could cause the observed rock sequence transition.

    The image shows a transition from a deeper water rock (shale, Rock A) to a
    shallower water rock (Rock B). This is a shallowing-upward sequence,
    also known as a marine regression, which is caused by a relative drop in sea level.

    Let's evaluate the options:
    i. Onset of glaciation: Locks up ocean water in ice sheets, causing a global sea-level drop. This leads to regression. (Plausible)
    ii. Coastal subsidence: The land sinks, causing a relative sea-level rise. This leads to transgression. (Not plausible)
    iii. Global warming: Melts ice and causes thermal expansion of water, leading to a global sea-level rise. This leads to transgression. (Not plausible)
    iv. Decreased Mid-Ocean Ridge Activity: Ridges cool and shrink, increasing the ocean basin volume and causing a global sea-level drop. This leads to regression. (Plausible)

    The plausible events are i and iv.
    """
    plausible_events = ['i', 'iv']
    # The required output format is a comma-separated list with no spaces.
    answer = ",".join(plausible_events)
    print(answer)

solve_geology_problem()