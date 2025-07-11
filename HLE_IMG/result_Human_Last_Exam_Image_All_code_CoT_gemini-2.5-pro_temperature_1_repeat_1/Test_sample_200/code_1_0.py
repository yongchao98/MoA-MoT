def solve_geology_question():
    """
    This function analyzes the geological events that could lead to the observed rock sequence.
    
    The transition from Rock A (likely a deep-water shale/mudstone) to Rock B (a shallower-water deposit)
    indicates a regression, which is a fall in relative sea level. We need to identify the events
    that cause sea-level to fall.

    Let's evaluate the options:
    i. Onset of glaciation: This locks up ocean water in ice sheets, causing global sea level to fall. This is a plausible cause for regression.
    ii. Coastal subsidence: This is the sinking of the land, which causes a relative sea-level rise (transgression), not a fall. This is not plausible.
    iii. Global warming: This melts ice and causes thermal expansion of water, leading to a global sea-level rise (transgression). This is not plausible.
    iv. Decreased Mid-Ocean Ridge Activity: Slower seafloor spreading causes the ridges to cool and contract, increasing the ocean basin's volume and causing global sea level to fall. This is a plausible cause for regression.

    Therefore, the plausible events are i and iv.
    """
    
    # List of plausible option labels
    plausible_options = ['i', 'iv']
    
    # Format the answer as a comma-separated string with no spaces
    answer = ",".join(plausible_options)
    
    print(answer)

solve_geology_question()