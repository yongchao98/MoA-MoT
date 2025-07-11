def calculate_vmt_reduction():
    """
    Calculates the percentage reduction in Car VMT from Baseline to Proposed scenarios.
    """
    # Common data for both scenarios
    share_of_total_trips = [0.48, 0.22, 0.14, 0.16]
    avg_trip_distance = [1, 4, 7.5, 20]
    distance_categories = ["0-3 miles", "3-5 miles", "5-10 miles", "10+ miles"]

    # --- Baseline Scenario Data (Table 1) ---
    baseline_car_1_occ = [0.25, 0.38, 0.40, 0.42]
    baseline_car_multi_occ = [0.28, 0.43, 0.45, 0.47]

    # --- Proposed Scenario Data (Table 2) ---
    proposed_car_1_occ = [0.01, 0.08, 0.14, 0.16]
    proposed_car_multi_occ = [0.02, 0.19, 0.27, 0.33]

    # --- VMT Calculation ---
    # VMT is calculated as a weighted average, assuming 100 total trips.
    # Formula: sum(share_of_trips * (car_1_occ_share + car_multi_occ_share) * avg_distance)
    
    baseline_vmt = 0
    proposed_vmt = 0
    
    print("Calculating Baseline Car VMT:")
    baseline_vmt_eq_parts = []
    for i in range(len(distance_categories)):
        total_car_share = baseline_car_1_occ[i] + baseline_car_multi_occ[i]
        vmt_part = share_of_total_trips[i] * total_car_share * avg_trip_distance[i]
        baseline_vmt += vmt_part
        baseline_vmt_eq_parts.append(f"({share_of_total_trips[i]} * ({baseline_car_1_occ[i]} + {baseline_car_multi_occ[i]}) * {avg_trip_distance[i]})")
    
    print("Baseline VMT Equation:")
    print(" + ".join(baseline_vmt_eq_parts))
    print(f"= {baseline_vmt:.4f}\n")

    print("Calculating Proposed Car VMT:")
    proposed_vmt_eq_parts = []
    for i in range(len(distance_categories)):
        total_car_share = proposed_car_1_occ[i] + proposed_car_multi_occ[i]
        vmt_part = share_of_total_trips[i] * total_car_share * avg_trip_distance[i]
        proposed_vmt += vmt_part
        proposed_vmt_eq_parts.append(f"({share_of_total_trips[i]} * ({proposed_car_1_occ[i]} + {proposed_car_multi_occ[i]}) * {avg_trip_distance[i]})")
        
    print("Proposed VMT Equation:")
    print(" + ".join(proposed_vmt_eq_parts))
    print(f"= {proposed_vmt:.4f}\n")

    # --- Percentage Reduction Calculation ---
    vmt_reduction = ((baseline_vmt - proposed_vmt) / baseline_vmt) * 100

    print("Calculating Percentage Reduction:")
    print(f"((Baseline VMT - Proposed VMT) / Baseline VMT) * 100")
    print(f"Percentage Reduction = (({baseline_vmt:.4f} - {proposed_vmt:.4f}) / {baseline_vmt:.4f}) * 100")
    print(f"Percentage Reduction = {vmt_reduction:.1f}%")

    return vmt_reduction

if __name__ == "__main__":
    final_reduction = calculate_vmt_reduction()
    print(f"<<<{final_reduction:.1f}>>>")