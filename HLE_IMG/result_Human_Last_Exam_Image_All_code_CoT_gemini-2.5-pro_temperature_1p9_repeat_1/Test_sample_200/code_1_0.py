def solve_geology_question():
    """
    Analyzes the geological events that could lead to the observed rock sequence.
    The transition from Rock A (deep water) to Rock B (shallower water) is a regression,
    which is caused by a relative fall in sea level.

    The function evaluates the given options:
    i. Onset of glaciation: Causes global sea level to fall. Plausible.
    ii. Coastal subsidence: Causes relative sea level to rise. Not plausible.
    iii. Global warming: Causes global sea level to rise. Not plausible.
    iv. Decreased Mid-Ocean Ridge Activity: Causes global sea level to fall. Plausible.
    """
    plausible_events = ['i', 'iv']
    answer = ",".join(plausible_events)
    print(answer)

solve_geology_question()
<<<i,iv>>>