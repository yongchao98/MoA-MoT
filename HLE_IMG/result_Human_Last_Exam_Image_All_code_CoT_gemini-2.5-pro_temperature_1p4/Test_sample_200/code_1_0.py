def solve_geology_problem():
    """
    This function analyzes the geological events that could cause a marine regression,
    as indicated by the transition from deep-water to shallow-water rock layers.
    """
    
    # The image shows a transition from a darker, likely deeper water rock (A)
    # to a lighter, likely shallower water rock (B). This is a regression,
    # which is caused by a fall in relative sea level.

    events = {
        'i': {'name': 'Onset of glaciation', 'effect': 'fall'},
        'ii': {'name': 'Coastal subsidence', 'effect': 'rise'},
        'iii': {'name': 'Global warming', 'effect': 'rise'},
        'iv': {'name': 'Decreased Mid-Ocean Ridge Activity', 'effect': 'fall'}
    }
    
    required_effect = 'fall'
    
    plausible_options = []
    for key, details in events.items():
        if details['effect'] == required_effect:
            plausible_options.append(key)
            
    # Sort the list for consistent output, e.g., ['i', 'iv']
    plausible_options.sort()
    
    # Format the answer as a comma-separated string with no spaces
    answer = ",".join(plausible_options)
    
    print(answer)

solve_geology_problem()
<<<i,iv>>>