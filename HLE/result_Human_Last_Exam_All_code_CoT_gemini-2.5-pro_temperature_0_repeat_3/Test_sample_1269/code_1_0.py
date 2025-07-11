def solve_hfr_mapping():
    """
    Calculates the transfer time for the 'azi' gene for different Hfr strains
    to find the scenario most consistent with early expression of 'azi'.
    """
    # Step 1: Define the E. coli gene map (positions in minutes)
    gene_map = {
        'thr': 0,
        'azi': 2,
        'ton': 3,
        'pro': 6,
        'lac': 8,
        'str': 73
    }

    # Define the scenarios from the answer choices
    scenarios = [
        {'label': 'A', 'origin': 'ton', 'direction': 'cw'},
        {'label': 'B', 'origin': 'lac', 'direction': 'ccw'},
        {'label': 'C', 'origin': 'pro', 'direction': 'cw'},
        {'label': 'D', 'origin': 'thr', 'direction': 'ccw'},
        {'label': 'E', 'origin': 'str', 'direction': 'cw'}
    ]

    target_gene = 'azi'
    target_pos = gene_map[target_gene]
    
    results = []

    print("Calculating the transfer time for the 'azi' gene in each scenario:\n")

    # Step 2: Calculate transfer time for each scenario
    for s in scenarios:
        origin_gene = s['origin']
        origin_pos = gene_map[origin_gene]
        direction = s['direction']
        time = 0
        
        if direction == 'cw': # Clockwise (increasing minutes)
            # Formula: (target_pos - origin_pos + 100) % 100
            time = (target_pos - origin_pos + 100) % 100
            print(f"Scenario {s['label']}: Origin near {origin_gene} ({origin_pos} min), Clockwise transfer.")
            print(f"  Equation: time = ({target_pos} - {origin_pos} + 100) % 100 = {time} minutes.")
        
        elif direction == 'ccw': # Counter-clockwise (decreasing minutes)
            # Formula: (origin_pos - target_pos + 100) % 100
            time = (origin_pos - target_pos + 100) % 100
            print(f"Scenario {s['label']}: Origin near {origin_gene} ({origin_pos} min), Counter-clockwise transfer.")
            print(f"  Equation: time = ({origin_pos} - {target_pos} + 100) % 100 = {time} minutes.")
        
        results.append({'label': s['label'], 'time': time})
        print("-" * 20)

    # Step 3: Identify the scenario with the minimum transfer time
    min_time_scenario = min(results, key=lambda x: x['time'])

    print("\nConclusion:")
    print(f"The scenario with the shortest transfer time for the 'azi' gene is Scenario {min_time_scenario['label']} with a time of {min_time_scenario['time']} minutes.")
    print("This is the most consistent with the experimental data showing early expression of the 'azi' gene.")

solve_hfr_mapping()
<<<B>>>