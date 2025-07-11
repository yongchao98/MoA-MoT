import pandas as pd

def calculate_vmt_reduction():
    """
    Calculates the percentage reduction in Car VMT based on the provided data tables.
    """
    # Data from the tables
    # Common data for both scenarios
    share_of_total_trips = [0.48, 0.22, 0.14, 0.16]
    average_trip_distance = [1, 4, 7.5, 20]
    distance_categories = ["0-3 miles", "3-5 miles", "5-10 miles", "10+ miles"]

    # Baseline scenario data
    baseline_car1_share = [0.25, 0.38, 0.40, 0.42]
    baseline_carpool_share = [0.28, 0.43, 0.45, 0.47]

    # Proposed scenario data
    proposed_car1_share = [0.01, 0.08, 0.14, 0.16]
    proposed_carpool_share = [0.02, 0.19, 0.27, 0.33]

    # --- Step 1: Calculate Baseline VMT ---
    baseline_vmt = 0
    print("Calculating Baseline Car VMT:")
    for i in range(len(distance_categories)):
        total_car_share = baseline_car1_share[i] + baseline_carpool_share[i]
        vmt_contribution = share_of_total_trips[i] * average_trip_distance[i] * total_car_share
        baseline_vmt += vmt_contribution
        print(f"  - Category '{distance_categories[i]}': {share_of_total_trips[i]:.2f} (trip share) * {average_trip_distance[i]:.1f} (miles) * ({baseline_car1_share[i]:.2f} + {baseline_carpool_share[i]:.2f}) (car share) = {vmt_contribution:.4f}")
    
    print(f"\nTotal Baseline Relative Car VMT = {baseline_vmt:.4f}\n")

    # --- Step 2: Calculate Proposed VMT ---
    proposed_vmt = 0
    print("Calculating Proposed Car VMT:")
    for i in range(len(distance_categories)):
        total_car_share = proposed_car1_share[i] + proposed_carpool_share[i]
        vmt_contribution = share_of_total_trips[i] * average_trip_distance[i] * total_car_share
        proposed_vmt += vmt_contribution
        print(f"  - Category '{distance_categories[i]}': {share_of_total_trips[i]:.2f} (trip share) * {average_trip_distance[i]:.1f} (miles) * ({proposed_car1_share[i]:.2f} + {proposed_carpool_share[i]:.2f}) (car share) = {vmt_contribution:.4f}")

    print(f"\nTotal Proposed Relative Car VMT = {proposed_vmt:.4f}\n")

    # --- Step 3: Calculate Percentage Reduction ---
    if baseline_vmt == 0:
        reduction_percentage = 0
    else:
        reduction_percentage = ((baseline_vmt - proposed_vmt) / baseline_vmt) * 100
    
    print("Calculating Percentage Reduction:")
    print(f"Formula: ((Baseline VMT - Proposed VMT) / Baseline VMT) * 100")
    print(f"Calculation: (({baseline_vmt:.4f} - {proposed_vmt:.4f}) / {baseline_vmt:.4f}) * 100 = {reduction_percentage:.1f}%")
    
    # Return the final numerical answer for the platform
    return reduction_percentage

# Run the calculation and store the final answer
final_answer = calculate_vmt_reduction()
# The final answer format is handled outside the function for clarity.
# <<<52.2>>>