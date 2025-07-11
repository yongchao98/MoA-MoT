def solve_geology_problem():
    """
    This function determines the plausible explanations for the observed rock strata transition.

    The transition from Rock A (dark, fine-grained, deep-water deposit) to Rock B 
    (lighter, coarser, shallower-water deposit) indicates a marine regression, or a fall in sea level.

    We evaluate the given options:
    i. Onset of glaciation: Causes a global sea-level fall. Plausible.
    ii. Coastal subsidence: Causes a relative sea-level rise (transgression). Not plausible.
    iii. Global warming: Causes a global sea-level rise (transgression). Not plausible.
    iv. Decreased Mid-Ocean Ridge Activity: Increases ocean basin volume, causing a global sea-level fall. Plausible.
    
    The plausible options are i and iv.
    """
    
    plausible_events = ['i', 'iv']
    
    # Format the answer as a comma-separated list with no spaces
    answer = ",".join(plausible_events)
    
    print(answer)

solve_geology_problem()