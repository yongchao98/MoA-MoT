def calculate_vmt_reduction():
    """
    Calculates the percentage reduction in Car VMT based on the provided tables.
    """
    # Data from the tables
    share_of_total_trips = [0.48, 0.22, 0.14, 0.16]
    avg_trip_distance = [1, 4, 7.5, 20]

    # Baseline mode shares for cars
    baseline_car_1_occ = [0.25, 0.38, 0.40, 0.42]
    baseline_car_multi_occ = [0.28, 0.43, 0.45, 0.47]

    # Proposed mode shares for cars
    proposed_car_1_occ = [0.01, 0.08, 0.14, 0.16]
    proposed_car_multi_occ = [0.02, 0.19, 0.27, 0.33]

    # --- Step 1: Calculate Baseline VMT ---
    baseline_vmt = 0
    print("Calculating Baseline VMT...")
    for i in range(len(share_of_total_trips)):
        # Total car share for the distance category
        total_car_share = baseline_car_1_occ[i] + baseline_car_multi_occ[i]
        # VMT for the distance category (Share of Trips * Car Mode Share * Avg Distance)
        vmt_component = share_of_total_trips[i] * total_car_share * avg_trip_distance[i]
        baseline_vmt += vmt_component
        print(f"  - Distance {i+1}: ({share_of_total_trips[i]} * ({baseline_car_1_occ[i]} + {baseline_car_multi_occ[i]}) * {avg_trip_distance[i]}) = {vmt_component:.4f}")
    print(f"Total Baseline VMT (relative value) = {baseline_vmt:.4f}\n")

    # --- Step 2: Calculate Proposed VMT ---
    proposed_vmt = 0
    print("Calculating Proposed VMT...")
    for i in range(len(share_of_total_trips)):
        # Total car share for the distance category
        total_car_share = proposed_car_1_occ[i] + proposed_car_multi_occ[i]
        # VMT for the distance category
        vmt_component = share_of_total_trips[i] * total_car_share * avg_trip_distance[i]
        proposed_vmt += vmt_component
        print(f"  - Distance {i+1}: ({share_of_total_trips[i]} * ({proposed_car_1_occ[i]} + {proposed_car_multi_occ[i]}) * {avg_trip_distance[i]}) = {vmt_component:.4f}")
    print(f"Total Proposed VMT (relative value) = {proposed_vmt:.4f}\n")

    # --- Step 3: Calculate Percentage Reduction ---
    reduction = baseline_vmt - proposed_vmt
    percentage_reduction = (reduction / baseline_vmt) * 100

    print("Calculating Percentage Reduction...")
    print(f"Reduction = Baseline VMT - Proposed VMT = {baseline_vmt:.4f} - {proposed_vmt:.4f} = {reduction:.4f}")
    print(f"Percentage Reduction = (Reduction / Baseline VMT) * 100")
    print(f"Percentage Reduction = ({reduction:.4f} / {baseline_vmt:.4f}) * 100 = {percentage_reduction:.1f}%")

    return percentage_reduction

# Run the calculation and print the final answer in the required format
final_answer = calculate_vmt_reduction()
# The final answer is wrapped in <<<...>>> as requested.
# print(f"\n<<<{final_answer:.1f}>>>")
# The problem asks to still output each number in the final equation. The lines above do that.
# The following is just the final numerical answer as requested.

print(f"\nBy changing the mode share, the Car VMT can be reduced by {final_answer:.1f}%.")
print(f"<<<{final_answer:.1f}>>>")