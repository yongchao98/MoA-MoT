import pandas as pd

def solve_vmt_reduction():
    """
    Calculates the percentage reduction in Car VMT based on the provided tables.
    """
    # Step 1: Define the data from the tables
    # Data is consistent across both tables unless specified
    distance_categories = ['0-3 miles', '3-5 miles', '5-10 miles', '10+ miles']
    share_of_total_trips = [0.48, 0.22, 0.14, 0.16]  # 48%, 22%, 14%, 16%
    average_trip_distance = [1, 4, 7.5, 20]

    # Baseline mode shares (Table 1)
    baseline_car_1_occ = [0.25, 0.38, 0.40, 0.42]
    baseline_car_gt2_occ = [0.28, 0.43, 0.45, 0.47]

    # Proposed mode shares (Table 2)
    proposed_car_1_occ = [0.01, 0.08, 0.14, 0.16]
    proposed_car_gt2_occ = [0.02, 0.19, 0.27, 0.33]

    # Combine car shares for total car mode share
    baseline_total_car_share = [a + b for a, b in zip(baseline_car_1_occ, baseline_car_gt2_occ)]
    proposed_total_car_share = [a + b for a, b in zip(proposed_car_1_occ, proposed_car_gt2_occ)]

    # Step 2 & 3: Calculate VMT Index for Baseline and Proposed
    baseline_vmt_terms = []
    proposed_vmt_terms = []

    for i in range(len(distance_categories)):
        baseline_vmt_terms.append(share_of_total_trips[i] * baseline_total_car_share[i] * average_trip_distance[i])
        proposed_vmt_terms.append(share_of_total_trips[i] * proposed_total_car_share[i] * average_trip_distance[i])
    
    baseline_vmt = sum(baseline_vmt_terms)
    proposed_vmt = sum(proposed_vmt_terms)

    # Print the calculation steps
    print("This calculation determines a VMT (Vehicle Miles Traveled) index, representing the average miles traveled by car per trip across the entire system.")
    print("\n--- Baseline Car VMT Index Calculation ---")
    baseline_eq = " + ".join([f"({sot:.2f} * {cs:.2f} * {d})" for sot, cs, d in zip(share_of_total_trips, baseline_total_car_share, average_trip_distance)])
    print(f"Baseline VMT Index = {baseline_eq}")
    print(f"                     = {' + '.join([f'{t:.4f}' for t in baseline_vmt_terms])}")
    print(f"                     = {baseline_vmt:.4f}")

    print("\n--- Proposed Car VMT Index Calculation ---")
    proposed_eq = " + ".join([f"({sot:.2f} * {cs:.2f} * {d})" for sot, cs, d in zip(share_of_total_trips, proposed_total_car_share, average_trip_distance)])
    print(f"Proposed VMT Index = {proposed_eq}")
    print(f"                   = {' + '.join([f'{t:.4f}' for t in proposed_vmt_terms])}")
    print(f"                   = {proposed_vmt:.4f}")

    # Step 4: Calculate the percentage reduction
    vmt_reduction = baseline_vmt - proposed_vmt
    percentage_reduction = (vmt_reduction / baseline_vmt) * 100

    print("\n--- Percentage Reduction Calculation ---")
    print(f"Reduction = (Baseline VMT - Proposed VMT) / Baseline VMT")
    print(f"          = (({baseline_vmt:.4f} - {proposed_vmt:.4f}) / {baseline_vmt:.4f}) * 100")
    print(f"          = ({vmt_reduction:.4f} / {baseline_vmt:.4f}) * 100")
    print(f"          = {percentage_reduction:.1f}%")
    
    return percentage_reduction

# Run the calculation and get the final answer
final_answer = solve_vmt_reduction()
print(f"\nBy changing the composition of mode share, the Car VMT can be reduced by {final_answer:.1f}%.")
print(f"<<<{final_answer:.1f}>>>")