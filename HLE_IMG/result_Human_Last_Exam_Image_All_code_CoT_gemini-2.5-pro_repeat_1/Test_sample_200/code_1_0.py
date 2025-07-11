def solve_geology_question():
    """
    This function determines the plausible explanations for the observed rock transition.

    The transition from the dark Rock A (likely deep-water shale) to the lighter Rock B 
    (likely shallower-water sediment) indicates a marine regression, which is caused by a 
    fall in relative sea level.

    Let's analyze the options:
    i. Onset of glaciation: Causes a global sea-level fall by locking water in ice sheets. This leads to regression. (Plausible)
    ii. Coastal subsidence: Causes a relative sea-level rise as the land sinks. This leads to transgression. (Not plausible)
    iii. Global warming: Causes a global sea-level rise from melting ice and thermal expansion. This leads to transgression. (Not plausible)
    iv. Decreased Mid-Ocean Ridge Activity: Causes a global sea-level fall as the ocean basin volume increases. This leads to regression. (Plausible)

    The plausible options are i and iv.
    """
    
    plausible_options = ['i', 'iv']
    
    # Format the answer as a comma-separated list with no spaces
    answer = ",".join(plausible_options)
    
    print(answer)

solve_geology_question()
<<<i,iv>>>