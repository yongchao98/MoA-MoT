def solve_geology_question():
    """
    Analyzes the geological scenario and determines the plausible causes.

    The transition from Rock A (dark, likely deep-water shale) to Rock B 
    (lighter, likely shallower-water sediment) indicates a fall in relative sea level, 
    a process called regression.

    We evaluate the options based on their effect on sea level:
    i. Onset of glaciation: Locks water in ice sheets, causing global sea level to fall (Regression). Plausible.
    ii. Coastal subsidence: Land sinks, causing relative sea level to rise (Transgression). Not plausible.
    iii. Global warming: Melts ice and thermal expansion, causing global sea level to rise (Transgression). Not plausible.
    iv. Decreased Mid-Ocean Ridge Activity: Ocean basins get larger, causing global sea level to fall (Regression). Plausible.
    """
    
    plausible_events = []

    # Event i: Onset of glaciation -> Regression
    plausible_events.append('i')

    # Event ii: Coastal subsidence -> Transgression (not correct)

    # Event iii: Global warming -> Transgression (not correct)

    # Event iv: Decreased Mid-Ocean Ridge Activity -> Regression
    plausible_events.append('iv')

    # Format the answer as a comma-separated list with no spaces
    answer = ",".join(plausible_events)
    print(answer)

solve_geology_question()