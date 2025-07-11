def calculate_vmt_reduction():
    """
    Calculates the percentage reduction in Car VMT based on the provided tables.
    """
    # Data from the tables
    share_of_total_trips = [0.48, 0.22, 0.14, 0.16]
    avg_trip_distance = [1, 4, 7.5, 20]
    distance_categories = ["0-3 miles", "3-5 miles", "5-10 miles", "10+ miles"]

    # Baseline car mode shares (1 occupant + carpool)
    baseline_car_1_occ = [0.25, 0.38, 0.40, 0.42]
    baseline_car_pool = [0.28, 0.43, 0.45, 0.47]
    baseline_total_car_share = [x + y for x, y in zip(baseline_car_1_occ, baseline_car_pool)]

    # Proposed car mode shares (1 occupant + carpool)
    proposed_car_1_occ = [0.01, 0.08, 0.14, 0.16]
    proposed_car_pool = [0.02, 0.19, 0.27, 0.33]
    proposed_total_car_share = [x + y for x, y in zip(proposed_car_1_occ, proposed_car_pool)]

    # --- Step 1: Calculate Baseline VMT Index ---
    print("--- Calculating Baseline Car VMT Index ---")
    baseline_vmt = 0
    for i in range(len(distance_categories)):
        vmt_component = share_of_total_trips[i] * baseline_total_car_share[i] * avg_trip_distance[i]
        print(f"Category '{distance_categories[i]}': {share_of_total_trips[i]} (trip share) * {baseline_total_car_share[i]:.2f} (car share) * {avg_trip_distance[i]} (miles) = {vmt_component:.4f}")
        baseline_vmt += vmt_component
    print(f"Total Baseline Car VMT Index = {baseline_vmt:.4f}\n")

    # --- Step 2: Calculate Proposed VMT Index ---
    print("--- Calculating Proposed Car VMT Index ---")
    proposed_vmt = 0
    for i in range(len(distance_categories)):
        vmt_component = share_of_total_trips[i] * proposed_total_car_share[i] * avg_trip_distance[i]
        print(f"Category '{distance_categories[i]}': {share_of_total_trips[i]} (trip share) * {proposed_total_car_share[i]:.2f} (car share) * {avg_trip_distance[i]} (miles) = {vmt_component:.4f}")
        proposed_vmt += vmt_component
    print(f"Total Proposed Car VMT Index = {proposed_vmt:.4f}\n")

    # --- Step 3: Calculate Percentage Reduction ---
    print("--- Calculating Percentage Reduction ---")
    reduction = baseline_vmt - proposed_vmt
    percentage_reduction = (reduction / baseline_vmt) * 100
    
    print(f"Percentage Reduction = (({baseline_vmt:.4f} - {proposed_vmt:.4f}) / {baseline_vmt:.4f}) * 100")
    print(f"Percentage Reduction = ({reduction:.4f} / {baseline_vmt:.4f}) * 100")
    print(f"Final Answer: The percentage of Car VMT that can be reduced is {percentage_reduction:.1f}%")

    return percentage_reduction

# Execute the calculation and print the result
final_answer = calculate_vmt_reduction()
# The final answer is wrapped in its own print statement for clarity, but the core logic is above.
# print(f"\nFinal calculated percentage reduction: {final_answer:.1f}%")

# Final answer in the required format
print(f"\n<<<{final_answer:.1f}>>>")