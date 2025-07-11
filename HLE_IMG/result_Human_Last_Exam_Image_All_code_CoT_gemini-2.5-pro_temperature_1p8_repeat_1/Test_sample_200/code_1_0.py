def solve_geology_question():
    """
    Analyzes the geological sequence and determines the plausible causes for the observed transition.

    The image shows a transition from a lower, darker rock layer (A) to an upper, lighter layer (B).
    Assuming standard deposition, A is older than B.
    - Rock A (dark shale/mudstone) suggests a low-energy, deep-water depositional environment.
    - Rock B (lighter, possibly siltstone/sandstone/limestone) suggests a higher-energy, shallower-water environment.
    - The transition from deep-water to shallow-water deposits is a regression, caused by a relative fall in sea level.

    We evaluate each option's effect on sea level:
    i. Onset of glaciation: Locks ocean water into ice sheets, causing global sea level to FALL. (Plausible)
    ii. Coastal subsidence: The land sinks, causing relative sea level to RISE. (Not plausible)
    iii. Global warming: Melts ice and causes thermal expansion of water, causing global sea level to RISE. (Not plausible)
    iv. Decreased Mid-Ocean Ridge Activity: Ocean ridges cool and shrink, increasing ocean basin volume, causing global sea level to FALL. (Plausible)

    Therefore, the plausible explanations are (i) and (iv).
    """
    plausible_events = ['i', 'iv']
    answer = ",".join(plausible_events)
    print(answer)

solve_geology_question()