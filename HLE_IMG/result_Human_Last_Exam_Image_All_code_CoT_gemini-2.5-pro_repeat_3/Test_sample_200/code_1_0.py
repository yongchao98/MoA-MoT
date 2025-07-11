def solve_geology_question():
    """
    This function determines the plausible explanations for the observed geological transition.
    The transition from Rock A (deeper water shale) to Rock B (shallower water sandstone/siltstone)
    represents a marine regression, or a fall in relative sea level.

    We evaluate the options:
    i. Onset of glaciation: Locks water in ice sheets, causing global sea level to fall. Plausible.
    ii. Coastal subsidence: Sinking of land, causes relative sea level to rise. Not plausible.
    iii. Global warming: Melts ice, causes global sea level to rise. Not plausible.
    iv. Decreased Mid-Ocean Ridge Activity: Cools and shrinks ridges, increasing ocean basin volume,
                                            causing global sea level to fall. Plausible.

    The correct options are those that cause a sea-level fall.
    """
    plausible_events = ['i', 'iv']
    answer = ",".join(plausible_events)
    print(answer)

solve_geology_question()