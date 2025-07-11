def solve_hfr_mapping():
    """
    Analyzes Hfr conjugation scenarios to find the one consistent
    with the early transfer of the 'azi' gene.
    """
    # Step 1: Define the E. coli gene map (positions in minutes)
    genes = {
        'thr': 0.0,
        'azi': 2.5,
        'ton': 3.2,
        'pro': 5.8,
        'lac': 8.0,
        'str': 73.0
    }
    
    # The experimental result indicates early transfer of 'azi'.
    # This means we are looking for the scenario with the minimum transfer time for 'azi'.
    print("The experimental data shows the 'azi' gene is transferred early.")
    print("This means the Hfr strain must have its origin of transfer (oriT) located")
    print("such that the distance to 'azi' in the direction of transfer is minimal.\n")
    print("Let's calculate the transfer time for 'azi' in each scenario:")
    print("-" * 60)

    # Step 2: Define the five scenarios from the answer choices
    scenarios = [
        {'name': 'A', 'origin_gene': 'ton', 'direction': 'Clockwise'},
        {'name': 'B', 'origin_gene': 'lac', 'direction': 'Counterclockwise'},
        {'name': 'C', 'origin_gene': 'pro', 'direction': 'Clockwise'},
        {'name': 'D', 'origin_gene': 'thr', 'direction': 'Counterclockwise'},
        {'name': 'E', 'origin_gene': 'str', 'direction': 'Clockwise'}
    ]

    azi_pos = genes['azi']
    results = []

    # Step 3: Loop through each scenario and calculate transfer time for 'azi'
    for sc in scenarios:
        origin_pos = genes[sc['origin_gene']]
        
        if sc['direction'] == 'Clockwise':
            # For clockwise transfer, time = (gene_pos - origin_pos) mod 100
            transfer_time = (azi_pos - origin_pos + 100) % 100
            equation_str = f"({azi_pos} - {origin_pos} + 100) % 100"
        else: # Counterclockwise
            # For counter-clockwise transfer, time = (origin_pos - gene_pos) mod 100
            transfer_time = (origin_pos - azi_pos + 100) % 100
            equation_str = f"({origin_pos} - {azi_pos} + 100) % 100"
        
        results.append({'name': sc['name'], 'time': transfer_time})
        
        print(f"Scenario {sc['name']}: Origin near {sc['origin_gene']} ({origin_pos}'), {sc['direction']}")
        print(f"  Calculation for 'azi' ({azi_pos}'): {equation_str} = {transfer_time:.1f} minutes")

    print("-" * 60)

    # Step 4: Find the best scenario
    best_scenario = min(results, key=lambda x: x['time'])
    
    print(f"\nConclusion:")
    print(f"Scenario {best_scenario['name']} has the shortest transfer time for the 'azi' gene ({best_scenario['time']:.1f} minutes).")
    print("While 'azi' is not the absolute first gene transferred in this scenario, it is part of the first group")
    print("of linked genes to be transferred, whereas in all other scenarios it is transferred very late.")
    print("This makes it the most consistent with the experimental observation.")


solve_hfr_mapping()
<<<B>>>