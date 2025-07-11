import sys

def solve_vmt_reduction():
    """
    Calculates the percentage reduction in Car VMT based on the provided tables.
    """
    # Step 1: Define data from the tables
    distance_categories = ['0-3 miles', '3-5 miles', '5-10 miles', '10+ miles']
    share_of_total_trips = [0.48, 0.22, 0.14, 0.16]
    avg_trip_distance = [1, 4, 7.5, 20]

    # Baseline scenario data (Table 1)
    baseline_car_1_occupant = [0.25, 0.38, 0.40, 0.42]
    baseline_car_multi_occupant = [0.28, 0.43, 0.45, 0.47]

    # Proposed scenario data (Table 2)
    proposed_car_1_occupant = [0.01, 0.08, 0.14, 0.16]
    proposed_car_multi_occupant = [0.02, 0.19, 0.27, 0.33]

    # Calculate total car mode shares for each scenario
    baseline_total_car_share = [s1 + s2 for s1, s2 in zip(baseline_car_1_occupant, baseline_car_multi_occupant)]
    proposed_total_car_share = [s1 + s2 for s1, s2 in zip(proposed_car_1_occupant, proposed_car_multi_occupant)]

    # --- Step 2: Calculate Baseline VMT Index ---
    baseline_vmt_components = []
    print("--- Calculating Baseline VMT Index ---")
    print("VMT Index = Sum of [ (Trip Share) * (Total Car Share) * (Avg Distance) ] for each category.\n")
    for i in range(len(distance_categories)):
        component = share_of_total_trips[i] * baseline_total_car_share[i] * avg_trip_distance[i]
        baseline_vmt_components.append(component)
        print(f"Category '{distance_categories[i]}':")
        print(f"  Car Share = {baseline_car_1_occupant[i]} + {baseline_car_multi_occupant[i]} = {baseline_total_car_share[i]:.2f}")
        print(f"  VMT Component = {share_of_total_trips[i]} * {baseline_total_car_share[i]:.2f} * {avg_trip_distance[i]} = {component:.4f}\n")
    
    baseline_vmt_index = sum(baseline_vmt_components)
    print(f"Total Baseline VMT Index = {' + '.join([f'{c:.4f}' for c in baseline_vmt_components])} = {baseline_vmt_index:.4f}\n")

    # --- Step 3: Calculate Proposed VMT Index ---
    proposed_vmt_components = []
    print("--- Calculating Proposed VMT Index ---")
    for i in range(len(distance_categories)):
        component = share_of_total_trips[i] * proposed_total_car_share[i] * avg_trip_distance[i]
        proposed_vmt_components.append(component)
        print(f"Category '{distance_categories[i]}':")
        print(f"  Car Share = {proposed_car_1_occupant[i]} + {proposed_car_multi_occupant[i]} = {proposed_total_car_share[i]:.2f}")
        print(f"  VMT Component = {share_of_total_trips[i]} * {proposed_total_car_share[i]:.2f} * {avg_trip_distance[i]} = {component:.4f}\n")
        
    proposed_vmt_index = sum(proposed_vmt_components)
    print(f"Total Proposed VMT Index = {' + '.join([f'{c:.4f}' for c in proposed_vmt_components])} = {proposed_vmt_index:.4f}\n")

    # --- Step 4: Calculate Percentage Reduction ---
    print("--- Calculating Percentage Reduction in VMT ---")
    vmt_reduction = baseline_vmt_index - proposed_vmt_index
    percentage_reduction = (vmt_reduction / baseline_vmt_index) * 100
    
    print("Final Equation: ((Baseline VMT - Proposed VMT) / Baseline VMT) * 100")
    print(f"Percentage Reduction = (({baseline_vmt_index:.4f} - {proposed_vmt_index:.4f}) / {baseline_vmt_index:.4f}) * 100")
    print(f"Percentage Reduction = ({vmt_reduction:.4f} / {baseline_vmt_index:.4f}) * 100 = {percentage_reduction:.1f}%")
    
    # Final answer format
    sys.stdout.write(f"\n<<<{percentage_reduction:.1f}>>>")

solve_vmt_reduction()