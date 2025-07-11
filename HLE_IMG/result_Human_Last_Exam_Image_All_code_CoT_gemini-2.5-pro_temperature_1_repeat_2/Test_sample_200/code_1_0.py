def solve_geology_puzzle():
    """
    This function determines the plausible explanations for the observed geological transition.
    
    The transition from Rock A (deep-water deposit) to Rock B (shallower-water deposit) represents a regression, which is caused by a fall in relative sea level.
    
    We evaluate the options:
    i. Onset of glaciation: Causes a global sea-level fall. Plausible.
    ii. Coastal subsidence: Causes a relative sea-level rise. Not plausible.
    iii. Global warming: Causes a global sea-level rise. Not plausible.
    iv. Decreased Mid-Ocean Ridge Activity: Increases ocean basin volume, causing a global sea-level fall. Plausible.
    
    The plausible explanations are i and iv.
    """
    
    # The list of plausible options identified
    plausible_options = ['i', 'iv']
    
    # Format the answer as a comma-separated string with no spaces
    answer = ",".join(plausible_options)
    
    print(answer)

solve_geology_puzzle()