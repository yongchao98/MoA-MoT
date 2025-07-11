import collections

def solve_geology_puzzle():
    """
    This function analyzes the geological events that could lead to the observed rock sequence.
    """
    
    # Step 1: Define the problem.
    # The transition from Rock A (deep-water deposit, likely shale) to Rock B (shallower-water deposit, likely sandstone)
    # represents a marine regression, which is a relative fall in sea level.
    # We need to find which events cause a relative sea-level fall.
    
    events = {
        'i': ('Onset of glaciation', 'fall'),
        'ii': ('Coastal subsidence', 'rise'),
        'iii': ('Global warming', 'rise'),
        'iv': ('Decreased Mid-Ocean Ridge Activity', 'fall')
    }
    
    print("Analysis:")
    print("The rock sequence shows a transition from a deeper water environment (Rock A) to a shallower water environment (Rock B).")
    print("This is called a marine regression, caused by a relative fall in sea level.")
    print("We must identify which of the following events cause a sea-level fall.")
    print("-" * 20)
    
    plausible_options = []
    
    # Step 2: Evaluate each event.
    for key, (description, effect) in events.items():
        is_plausible = (effect == 'fall')
        plausibility_text = "plausible" if is_plausible else "not plausible"
        causation_text = "causes a sea-level fall" if effect == 'fall' else "causes a sea-level rise"
        
        print(f"Event {key}: {description}")
        print(f"This event {causation_text}, which is consistent with a regression.")
        print(f"Therefore, this option is {plausibility_text}.")
        print("-" * 20)
        
        if is_plausible:
            plausible_options.append(key)
            
    # Step 3: Format and print the final answer.
    final_answer = ",".join(plausible_options)
    print("Final Answer:")
    print(final_answer)

solve_geology_puzzle()