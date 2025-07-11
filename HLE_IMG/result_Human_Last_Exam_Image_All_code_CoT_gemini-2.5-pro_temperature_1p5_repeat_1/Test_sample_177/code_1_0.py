def calculate_vmt_reduction():
    """
    Calculates the percentage reduction in Car VMT based on the provided tables.
    """
    # Data from the tables
    # Common data for both scenarios
    share_of_total_trips = [0.48, 0.22, 0.14, 0.16]
    avg_trip_distance = [1, 4, 7.5, 20]
    distance_categories = ["0-3 miles", "3-5 miles", "5-10 miles", "10+ miles"]

    # Table 1: Baseline mode shares for cars
    baseline_car_1_occ = [0.25, 0.38, 0.40, 0.42]
    baseline_car_multi_occ = [0.28, 0.43, 0.45, 0.47]

    # Table 2: Proposed mode shares for cars
    proposed_car_1_occ = [0.01, 0.08, 0.14, 0.16]
    proposed_car_multi_occ = [0.02, 0.19, 0.27, 0.33]

    def calculate_total_vmt(car_1_occ_shares, car_multi_occ_shares):
        """Helper function to calculate total VMT for a given scenario."""
        total_vmt = 0
        for i in range(len(distance_categories)):
            # Total share of trips made by car in this distance category
            total_car_share = car_1_occ_shares[i] + car_multi_occ_shares[i]
            # VMT for this category is (share of all trips * total car share * avg distance)
            vmt_for_category = share_of_total_trips[i] * total_car_share * avg_trip_distance[i]
            total_vmt += vmt_for_category
        return total_vmt

    # 1. Calculate Baseline VMT
    baseline_vmt = calculate_total_vmt(baseline_car_1_occ, baseline_car_multi_occ)
    print("Calculating Baseline Car VMT...")
    # Multiplying by 100 for a more intuitive number (VMT per 100 total trips)
    print(f"Total Baseline Car VMT: {baseline_vmt * 100:.2f}")
    print("-" * 30)

    # 2. Calculate Proposed VMT
    proposed_vmt = calculate_total_vmt(proposed_car_1_occ, proposed_car_multi_occ)
    print("Calculating Proposed Car VMT...")
    print(f"Total Proposed Car VMT: {proposed_vmt * 100:.2f}")
    print("-" * 30)
    
    # 3. Calculate Percentage Reduction
    vmt_reduction = baseline_vmt - proposed_vmt
    percentage_reduction = (vmt_reduction / baseline_vmt) * 100

    print("Calculating Percentage Reduction in Car VMT:")
    print(f"Equation: (({baseline_vmt:.4f} - {proposed_vmt:.4f}) / {baseline_vmt:.4f}) * 100")
    print(f"Result: {percentage_reduction:.1f}%")

    return percentage_reduction

# Run the calculation and print the final answer in the required format
final_answer = calculate_vmt_reduction()
print(f"<<<{final_answer:.1f}>>>")
