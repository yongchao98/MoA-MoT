def find_plausible_causes():
    """
    This function determines the plausible geological events that could cause the
    transition observed from Rock A to Rock B, which signifies a marine regression
    (sea level fall).
    """

    # The transition from deep-water rock (A) to shallow-water rock (B) indicates a fall in sea level.
    # We will identify which options cause a sea level fall.
    
    options = {
        'i': 'Onset of glaciation',
        'ii': 'Coastal subsidence',
        'iii': 'Global warming',
        'iv': 'Decreased Mid-Ocean Ridge Activity'
    }

    # Analyzing the effect of each option on sea level.
    # 'fall' means it's a plausible cause for the observed regression.
    # 'rise' means it's not a plausible cause.
    effects = {
        'i': 'fall',
        'ii': 'rise',
        'iii': 'rise',
        'iv': 'fall'
    }

    print("Analysis of Geological Events:")
    print("The transition from Rock A (deep water) to Rock B (shallow water) indicates a sea level fall (regression).")
    print("-" * 30)

    plausible_events = []
    for key, description in options.items():
        effect = effects[key]
        if effect == 'fall':
            plausible_events.append(key)
            print(f"Event '{key}. {description}' causes a sea level FALL. -> Plausible")
        else:
            print(f"Event '{key}. {description}' causes a sea level RISE. -> Not Plausible")

    # Format the final answer as a comma-separated list with no spaces.
    final_answer = ",".join(plausible_events)
    
    print("-" * 30)
    print("The plausible explanations are i and iv.")
    print("Final Answer formatted as a comma-separated list:")
    print(final_answer)

find_plausible_causes()