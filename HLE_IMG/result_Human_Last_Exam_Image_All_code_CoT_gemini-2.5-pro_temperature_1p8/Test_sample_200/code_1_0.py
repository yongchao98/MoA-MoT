def solve_geology_question():
    """
    Analyzes geological events and their effect on sea level to explain a rock sequence.
    """
    # The rock sequence shows a transition from deep-water deposits (A) to
    # shallower-water deposits (B), indicating a fall in sea level (regression).
    # We need to find the events that cause a sea-level fall.
    
    events = {
        'i': "Onset of glaciation",
        'ii': "Coastal subsidence",
        'iii': "Global warming",
        'iv': "Decreased Mid-Ocean Ridge Activity"
    }
    
    # Effect of each event on sea level
    effects = {
        'i': "fall",  # Glaciation locks water in ice sheets
        'ii': "rise",  # Land sinks, relative sea level rises
        'iii': "rise",  # Ice melts, water expands
        'iv': "fall"   # Ocean basin volume increases
    }
    
    required_effect = "fall"
    
    correct_options = []
    for option, effect in effects.items():
        if effect == required_effect:
            correct_options.append(option)
            
    # Sort to maintain a consistent order (i, iv)
    correct_options.sort()
    
    # Format the answer as a comma-separated list with no spaces
    answer = ",".join(correct_options)
    
    print(answer)

solve_geology_question()