def calculate_vmt_reduction():
    """
    Calculates the percentage reduction in Car VMT based on the provided tables.
    """
    # Data from the tables
    share_of_trips = [0.48, 0.22, 0.14, 0.16]
    avg_distances = [1, 4, 7.5, 20]
    distance_categories = ["0-3 miles", "3-5 miles", "5-10 miles", "10+ miles"]

    # Baseline car mode shares
    baseline_car_1_occ = [0.25, 0.38, 0.40, 0.42]
    baseline_car_multi_occ = [0.28, 0.43, 0.45, 0.47]

    # Proposed car mode shares
    proposed_car_1_occ = [0.01, 0.08, 0.14, 0.16]
    proposed_car_multi_occ = [0.02, 0.19, 0.27, 0.33]

    # --- Step 1: Calculate Baseline Car VMT ---
    print("Step 1: Calculate Baseline Car VMT Index")
    print("VMT Index = Sum of (Share of Trips * Total Car Share * Average Distance) for each category\n")
    
    baseline_total_car_share = [(c1 + c2) for c1, c2 in zip(baseline_car_1_occ, baseline_car_multi_occ)]
    baseline_vmt_components = []
    
    print("Baseline Calculation:")
    for i in range(len(distance_categories)):
        component = share_of_trips[i] * baseline_total_car_share[i] * avg_distances[i]
        baseline_vmt_components.append(component)
        print(f"- {distance_categories[i]:<10}: {share_of_trips[i]} (trip share) * ({baseline_car_1_occ[i]} + {baseline_car_multi_occ[i]}) (car share) * {avg_distances[i]} (miles) = {component:.4f}")

    baseline_vmt = sum(baseline_vmt_components)
    print(f"\nTotal Baseline VMT Index = {' + '.join([f'{c:.4f}' for c in baseline_vmt_components])} = {baseline_vmt:.4f}\n")

    # --- Step 2: Calculate Proposed Car VMT ---
    print("Step 2: Calculate Proposed Car VMT Index")
    
    proposed_total_car_share = [(c1 + c2) for c1, c2 in zip(proposed_car_1_occ, proposed_car_multi_occ)]
    proposed_vmt_components = []

    print("Proposed Calculation:")
    for i in range(len(distance_categories)):
        component = share_of_trips[i] * proposed_total_car_share[i] * avg_distances[i]
        proposed_vmt_components.append(component)
        print(f"- {distance_categories[i]:<10}: {share_of_trips[i]} (trip share) * ({proposed_car_1_occ[i]} + {proposed_car_multi_occ[i]}) (car share) * {avg_distances[i]} (miles) = {component:.4f}")

    proposed_vmt = sum(proposed_vmt_components)
    print(f"\nTotal Proposed VMT Index = {' + '.join([f'{c:.4f}' for c in proposed_vmt_components])} = {proposed_vmt:.4f}\n")

    # --- Step 3: Calculate Percentage Reduction ---
    print("Step 3: Calculate Percentage Reduction")
    vmt_reduction = baseline_vmt - proposed_vmt
    percentage_reduction = (vmt_reduction / baseline_vmt) * 100
    
    print(f"VMT Reduction = Baseline VMT - Proposed VMT")
    print(f"VMT Reduction = {baseline_vmt:.4f} - {proposed_vmt:.4f} = {vmt_reduction:.4f}")
    
    print(f"\nPercentage Reduction = (VMT Reduction / Baseline VMT) * 100")
    print(f"Percentage Reduction = ({vmt_reduction:.4f} / {baseline_vmt:.4f}) * 100 = {percentage_reduction:.1f}%")
    
    return percentage_reduction

if __name__ == '__main__':
    result = calculate_vmt_reduction()
    print(f"\n<<<The percentage of Car VMT that can be reduced is {result:.1f}%>>>")
    print(f"<<<{result:.1f}>>>")