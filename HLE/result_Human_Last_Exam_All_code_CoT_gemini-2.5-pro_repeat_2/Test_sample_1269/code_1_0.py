def solve_hfr_mapping():
    """
    Calculates the transfer time for the azis gene for various Hfr strains
    to determine which is most consistent with early expression.
    """
    # 1. Establish a standard E. coli gene map (positions in minutes)
    gene_map = {
        'thr': 0,
        'azi': 2,
        'ton': 2.5,
        'pro': 6,
        'lac': 8,
        'str': 73,
    }
    target_gene = 'azi'
    chromosome_size = 100

    # 2. Define the Hfr strains from the answer choices
    options = [
        {'id': 'A', 'origin_gene': 'ton', 'direction': 'Clockwise'},
        {'id': 'B', 'origin_gene': 'lac', 'direction': 'Counterclockwise'},
        {'id': 'C', 'origin_gene': 'pro', 'direction': 'Clockwise'},
        {'id': 'D', 'origin_gene': 'thr', 'direction': 'Counterclockwise'},
        {'id': 'E', 'origin_gene': 'str', 'direction': 'Clockwise'},
    ]

    print("Analyzing Hfr strains to find which transfers the 'azi' gene earliest.\n")
    print(f"Gene positions (minutes): {gene_map}")
    print(f"Target gene: '{target_gene}' at {gene_map[target_gene]} minutes.")
    print("-" * 60)

    results = {}

    # 3. Calculate transfer time for each option
    for option in options:
        origin_pos = gene_map[option['origin_gene']]
        target_pos = gene_map[target_gene]
        direction = option['direction']
        time = 0

        if direction == 'Clockwise':
            if target_pos >= origin_pos:
                time = target_pos - origin_pos
            else: # Wraps around the circle
                time = (chromosome_size - origin_pos) + target_pos
        
        elif direction == 'Counterclockwise':
            if target_pos <= origin_pos:
                time = origin_pos - target_pos
            else: # Wraps around the circle
                time = origin_pos + (chromosome_size - target_pos)

        results[option['id']] = time
        print(f"Option {option['id']}: Origin near '{option['origin_gene']}' ({origin_pos} min), Direction: {direction}")
        print(f"   Calculation: The time to transfer '{target_gene}' ({target_pos} min) is {time:.1f} minutes.")
        
    # 4. Find the best option (minimum transfer time)
    best_option_id = min(results, key=results.get)
    
    print("-" * 60)
    print(f"Conclusion: Option {best_option_id} results in the shortest transfer time for the 'azi' gene ({results[best_option_id]:.1f} minutes).")
    print("This is the most consistent scenario with the experimental data.")

solve_hfr_mapping()
<<<B>>>