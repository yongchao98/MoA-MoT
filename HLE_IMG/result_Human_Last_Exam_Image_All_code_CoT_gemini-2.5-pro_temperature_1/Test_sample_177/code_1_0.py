import pandas as pd

def calculate_vmt_reduction():
    """
    Calculates the percentage reduction in Car VMT based on the provided tables.
    """
    # Common data for both scenarios
    share_of_total_trips = [0.48, 0.22, 0.14, 0.16]
    avg_trip_distance = [1, 4, 7.5, 20]
    distance_bins = ["0-3 miles", "3-5 miles", "5-10 miles", "10+ miles"]

    # Baseline data (Table 1)
    baseline_car_1_occupant = [0.25, 0.38, 0.40, 0.42]
    baseline_car_multi_occupant = [0.28, 0.43, 0.45, 0.47]

    # Proposed data (Table 2)
    proposed_car_1_occupant = [0.01, 0.08, 0.14, 0.16]
    proposed_car_multi_occupant = [0.02, 0.19, 0.27, 0.33]

    # --- Step 1: Calculate Baseline Car VMT ---
    print("--- Calculating Baseline Car VMT ---")
    baseline_vmt_components = []
    for i in range(len(distance_bins)):
        total_car_share = baseline_car_1_occupant[i] + baseline_car_multi_occupant[i]
        vmt_component = share_of_total_trips[i] * total_car_share * avg_trip_distance[i]
        baseline_vmt_components.append(vmt_component)
    
    baseline_total_vmt = sum(baseline_vmt_components)
    print(f"Total Baseline Car VMT Index = {baseline_total_vmt:.4f}\n")

    # --- Step 2: Calculate Proposed Car VMT ---
    print("--- Calculating Proposed Car VMT ---")
    proposed_vmt_components = []
    for i in range(len(distance_bins)):
        total_car_share = proposed_car_1_occupant[i] + proposed_car_multi_occupant[i]
        vmt_component = share_of_total_trips[i] * total_car_share * avg_trip_distance[i]
        proposed_vmt_components.append(vmt_component)

    proposed_total_vmt = sum(proposed_vmt_components)
    print(f"Total Proposed Car VMT Index = {proposed_total_vmt:.4f}\n")

    # --- Step 3: Calculate Percentage Reduction ---
    print("--- Calculating Percentage Reduction ---")
    vmt_reduction = baseline_total_vmt - proposed_total_vmt
    percentage_reduction = (vmt_reduction / baseline_total_vmt) * 100

    print("Final Equation:")
    print(f"(({baseline_total_vmt:.4f} - {proposed_total_vmt:.4f}) / {baseline_total_vmt:.4f}) * 100 = {percentage_reduction:.1f}%")
    
    return percentage_reduction

# Run the calculation and store the final answer
final_answer = calculate_vmt_reduction()