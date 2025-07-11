import sys

def solve_hfr_mapping():
    """
    Calculates the transfer time for the azis gene in different Hfr scenarios
    to determine which is most consistent with early expression.
    """

    # Step 1: Establish the standard E. coli gene map positions in minutes.
    gene_map = {
        'thr': 0,
        'azi': 2,
        'ton': 3,
        'pro': 6,
        'lac': 8,
        'str': 73
    }
    
    total_map_minutes = 100
    target_gene = 'azi'
    target_pos = gene_map[target_gene]

    print("Analyzing Hfr transfer scenarios based on the E. coli gene map (100 minutes total).\n")
    print(f"Goal: Find the scenario where the '{target_gene}' gene (at {target_pos} min) is transferred the fastest.\n")

    results = {}

    # --- Scenario A: Clockwise, origin near ton ---
    origin_gene_a = 'ton'
    origin_pos_a = gene_map[origin_gene_a]
    # Clockwise transfer from 3 to 2 means wrapping around the circle.
    time_a = (total_map_minutes - origin_pos_a) + target_pos
    results['A'] = time_a
    print(f"A. Clockwise, origin near {origin_gene_a} (at {origin_pos_a} min):")
    print(f"   Transfer path is clockwise, wrapping around the map.")
    print(f"   Calculation: ({total_map_minutes} - {origin_pos_a}) + {target_pos} = {time_a} minutes.\n")

    # --- Scenario B: Counterclockwise, origin near lac ---
    origin_gene_b = 'lac'
    origin_pos_b = gene_map[origin_gene_b]
    # Counterclockwise transfer from 8 to 2 is a direct path.
    time_b = origin_pos_b - target_pos
    results['B'] = time_b
    print(f"B. Counterclockwise, origin near {origin_gene_b} (at {origin_pos_b} min):")
    print(f"   Transfer path is counterclockwise.")
    print(f"   Calculation: {origin_pos_b} - {target_pos} = {time_b} minutes.\n")

    # --- Scenario C: Clockwise, origin near pro ---
    origin_gene_c = 'pro'
    origin_pos_c = gene_map[origin_gene_c]
    # Clockwise transfer from 6 to 2 means wrapping around the circle.
    time_c = (total_map_minutes - origin_pos_c) + target_pos
    results['C'] = time_c
    print(f"C. Clockwise, origin near {origin_gene_c} (at {origin_pos_c} min):")
    print(f"   Transfer path is clockwise, wrapping around the map.")
    print(f"   Calculation: ({total_map_minutes} - {origin_pos_c}) + {target_pos} = {time_c} minutes.\n")

    # --- Scenario D: Counterclockwise, origin near thr ---
    origin_gene_d = 'thr'
    origin_pos_d = gene_map[origin_gene_d]
    # Counterclockwise transfer from 0 to 2 means wrapping around the circle.
    time_d = (origin_pos_d - target_pos) + total_map_minutes
    results['D'] = time_d
    print(f"D. Counterclockwise, origin near {origin_gene_d} (at {origin_pos_d} min):")
    print(f"   Transfer path is counterclockwise, wrapping around the map.")
    print(f"   Calculation: ({origin_pos_d} - {target_pos}) + {total_map_minutes} = {time_d} minutes.\n")

    # --- Scenario E: Clockwise, origin near str ---
    origin_gene_e = 'str'
    origin_pos_e = gene_map[origin_gene_e]
    # Clockwise transfer from 73 to 2 means wrapping around the circle.
    time_e = (total_map_minutes - origin_pos_e) + target_pos
    results['E'] = time_e
    print(f"E. Clockwise, origin near {origin_gene_e} (at {origin_pos_e} min):")
    print(f"   Transfer path is clockwise, wrapping around the map.")
    print(f"   Calculation: ({total_map_minutes} - {origin_pos_e}) + {target_pos} = {time_e} minutes.\n")

    # Find the best option (minimum time)
    best_option = min(results, key=results.get)
    min_time = results[best_option]

    print("--- Conclusion ---")
    print(f"The shortest transfer time for the '{target_gene}' gene is {min_time} minutes, which occurs in scenario {best_option}.")
    print("This scenario is the most consistent with the experimental observation of the gene being expressed early.")

solve_hfr_mapping()