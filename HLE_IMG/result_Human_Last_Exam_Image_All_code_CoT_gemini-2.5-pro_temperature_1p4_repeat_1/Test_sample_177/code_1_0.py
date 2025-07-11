def calculate_vmt_reduction():
    """
    Calculates the percentage reduction in Car VMT based on the provided tables.
    """
    # Data from the tables
    # Common data for both scenarios
    share_of_total_trips = [0.48, 0.22, 0.14, 0.16]
    average_trip_distance = [1, 4, 7.5, 20]

    # Baseline scenario data (Table 1)
    baseline_car_1_occupant_share = [0.25, 0.38, 0.40, 0.42]
    baseline_car_gt2_occupants_share = [0.28, 0.43, 0.45, 0.47]

    # Proposed scenario data (Table 2)
    proposed_car_1_occupant_share = [0.01, 0.08, 0.14, 0.16]
    proposed_car_gt2_occupants_share = [0.02, 0.19, 0.27, 0.33]

    # Step 1: Calculate Baseline VMT
    baseline_vmt = 0
    for i in range(len(share_of_total_trips)):
        total_car_share = baseline_car_1_occupant_share[i] + baseline_car_gt2_occupants_share[i]
        baseline_vmt += share_of_total_trips[i] * total_car_share * average_trip_distance[i]

    # Step 2: Calculate Proposed VMT
    proposed_vmt = 0
    for i in range(len(share_of_total_trips)):
        total_car_share = proposed_car_1_occupant_share[i] + proposed_car_gt2_occupants_share[i]
        proposed_vmt += share_of_total_trips[i] * total_car_share * average_trip_distance[i]

    # Step 3: Calculate Percentage Reduction
    vmt_reduction = baseline_vmt - proposed_vmt
    percentage_reduction = (vmt_reduction / baseline_vmt) * 100

    print("Step 1: Calculate Baseline Car VMT")
    print(f"Total Baseline Car VMT (weighted average) = {baseline_vmt:.4f} miles")
    print("\nStep 2: Calculate Proposed Car VMT")
    print(f"Total Proposed Car VMT (weighted average) = {proposed_vmt:.4f} miles")
    print("\nStep 3: Calculate Percentage Reduction")
    print("Percentage Reduction = ((Baseline VMT - Proposed VMT) / Baseline VMT) * 100")
    print(f"Percentage Reduction = (({baseline_vmt:.4f} - {proposed_vmt:.4f}) / {baseline_vmt:.4f}) * 100")
    print(f"\nThe percentage of Car VMT that can be reduced is: {percentage_reduction:.1f}%")

calculate_vmt_reduction()