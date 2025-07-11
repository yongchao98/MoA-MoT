def solve_conjugation_problem():
    """
    Calculates the transfer time for the 'azis' gene in different Hfr conjugation scenarios
    to determine which is most consistent with early expression.
    """
    # Step 1: Define the gene positions on the 100-minute E. coli chromosome map.
    genes = {
        'thr': 0,
        'leu': 2,
        'azis': 2.5,
        'ton': 3,
        'pro': 6,
        'lac': 8,
        'gal': 17,
        'his': 44,
        'str': 73
    }
    
    target_gene = 'azis'
    target_pos = genes[target_gene]
    
    # Step 2: Define the scenarios from the answer choices.
    scenarios = {
        'A': {'direction': 'Clockwise', 'origin_gene': 'ton'},
        'B': {'direction': 'Counterclockwise', 'origin_gene': 'lac'},
        'C': {'direction': 'Clockwise', 'origin_gene': 'pro'},
        'D': {'direction': 'Counterclockwise', 'origin_gene': 'thr'},
        'E': {'direction': 'Clockwise', 'origin_gene': 'str'}
    }
    
    print("### Analysis of Hfr Conjugation Scenarios ###")
    print(f"Goal: Find the scenario where the '{target_gene}' gene (at {target_pos} min) is transferred earliest.\n")

    results = {}

    # Step 3 & 4: Calculate the transfer time for 'azis' in each scenario.
    for key, params in scenarios.items():
        direction = params['direction']
        origin_gene = params['origin_gene']
        origin_pos = genes[origin_gene]
        
        print(f"Scenario {key}: {direction} transfer from an origin near '{origin_gene}' (at {origin_pos} min)")
        
        time = 0
        equation_str = ""
        
        if direction == 'Clockwise':
            if target_pos >= origin_pos:
                time = target_pos - origin_pos
                equation_str = f"Time = {target_pos} - {origin_pos}"
            else: # Wraps around the circle
                time = (100 - origin_pos) + target_pos
                equation_str = f"Time = (100 - {origin_pos}) + {target_pos}"
        
        elif direction == 'Counterclockwise':
            if target_pos <= origin_pos:
                time = origin_pos - target_pos
                equation_str = f"Time = {origin_pos} - {target_pos}"
            else: # Wraps around the circle
                time = origin_pos + (100 - target_pos)
                equation_str = f"Time = {origin_pos} + (100 - {target_pos})"
        
        results[key] = time
        print(f"  Calculation: {equation_str} = {time:.1f} minutes")
        print("-" * 20)

    # Step 5: Compare results and find the best option.
    best_option = min(results, key=results.get)
    min_time = results[best_option]
    
    print("\n### Conclusion ###")
    print(f"Comparing the calculated transfer times, the minimum time is {min_time:.1f} minutes, which occurs in Scenario {best_option}.")
    print("This scenario is the most consistent with the experimental result that the 'azis' gene is expressed before others.")

solve_conjugation_problem()
<<<B>>>