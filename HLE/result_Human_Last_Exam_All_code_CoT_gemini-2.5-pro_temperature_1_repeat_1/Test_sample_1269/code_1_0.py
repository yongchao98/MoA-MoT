def solve_hfr_mapping():
    """
    Calculates the transfer time for the 'azis' gene for different Hfr strains
    to determine which is most consistent with early expression.
    """
    # Step 1: Define gene locations on the 100-minute E. coli map.
    gene_map = {
        'thr': 0,
        'azis': 2,
        'pro': 6,
        'lac': 8,
        'ton': 10,
        'str': 73
    }

    target_gene_name = 'azis'
    target_gene_pos = gene_map[target_gene_name]
    chromosome_length = 100

    # Step 2: Define the scenarios from the answer choices.
    options = {
        'A': {'origin_gene': 'ton', 'direction': 'Clockwise'},
        'B': {'origin_gene': 'lac', 'direction': 'Counterclockwise'},
        'C': {'origin_gene': 'pro', 'direction': 'Clockwise'},
        'D': {'origin_gene': 'thr', 'direction': 'Counterclockwise'},
        'E': {'origin_gene': 'str', 'direction': 'Clockwise'}
    }

    print("Plan: Calculate the transfer time for the 'azis' gene (at 2 min) for each scenario.\n")
    
    results = {}

    # Step 3: Calculate transfer time for each option.
    for label, details in options.items():
        origin_gene = details['origin_gene']
        origin_pos = gene_map[origin_gene]
        direction = details['direction']
        time = 0
        
        equation = f"{label}. {direction} from {origin_gene} ({origin_pos} min): "
        equation += f"Transfer time to {target_gene_name} ({target_gene_pos} min) = "

        if direction == 'Clockwise':
            if target_gene_pos > origin_pos:
                time = target_gene_pos - origin_pos
                equation += f"{target_gene_pos} - {origin_pos} = {time}"
            else: # Wraps around
                time = (chromosome_length - origin_pos) + target_gene_pos
                equation += f"({chromosome_length} - {origin_pos}) + {target_gene_pos} = {time}"
        else: # Counterclockwise
            if origin_pos > target_gene_pos:
                time = origin_pos - target_gene_pos
                equation += f"{origin_pos} - {target_gene_pos} = {time}"
            else: # Wraps around
                time = (chromosome_length - target_gene_pos) + origin_pos
                equation += f"({chromosome_length} - {target_gene_pos}) + {origin_pos} = {time}"
        
        print(f"{equation} minutes")
        results[label] = time

    # Step 4: Find the best option (minimum transfer time).
    best_option = min(results, key=results.get)
    
    print("\nConclusion:")
    print(f"The scenario with the shortest transfer time for the 'azis' gene is Option {best_option}.")
    print("This is the most consistent with the experimental observation of prolonged expression of 'azis' before other genes.")

solve_hfr_mapping()
<<<B>>>