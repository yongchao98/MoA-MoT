import numpy as np

def calculate_vmt_reduction():
    """
    Calculates the percentage reduction in Car VMT from a Baseline to a Proposed scenario.
    """
    # Common data for both scenarios from the tables
    share_of_total_trips = np.array([0.48, 0.22, 0.14, 0.16])
    avg_trip_distance = np.array([1, 4, 7.5, 20])
    distance_labels = ["0-3 miles", "3-5 miles", "5-10 miles", "10+ miles"]

    # --- Baseline Scenario Data ---
    # Mode share for 'car, 1 occupant'
    baseline_car_1_occ = np.array([0.25, 0.38, 0.40, 0.42])
    # Mode share for 'car, with more than 2 occupants'
    baseline_car_multi_occ = np.array([0.28, 0.43, 0.45, 0.47])
    # Total car share is the sum of the two car modes
    baseline_total_car_share = baseline_car_1_occ + baseline_car_multi_occ

    # --- Proposed Scenario Data ---
    # Mode share for 'car, 1 occupant'
    proposed_car_1_occ = np.array([0.01, 0.08, 0.14, 0.16])
    # Mode share for 'car, with more than 2 occupants'
    proposed_car_multi_occ = np.array([0.02, 0.19, 0.27, 0.33])
    # Total car share is the sum of the two car modes
    proposed_total_car_share = proposed_car_1_occ + proposed_car_multi_occ

    # --- Calculations ---
    # VMT is proportional to: Share of Trips * Car Mode Share * Average Distance
    baseline_vmt_by_segment = share_of_total_trips * baseline_total_car_share * avg_trip_distance
    proposed_vmt_by_segment = share_of_total_trips * proposed_total_car_share * avg_trip_distance

    total_baseline_vmt = np.sum(baseline_vmt_by_segment)
    total_proposed_vmt = np.sum(proposed_vmt_by_segment)

    # Percentage reduction calculation
    vmt_reduction_percentage = ((total_baseline_vmt - total_proposed_vmt) / total_baseline_vmt) * 100

    # --- Output Results ---
    print("--- Calculating VMT Reduction ---")
    print("\nStep 1: Calculate Total Relative VMT for Baseline Scenario")
    vmt_baseline_eq = " + ".join([f"({s:.2f}*{c:.2f}*{d})" for s, c, d in zip(share_of_total_trips, baseline_total_car_share, avg_trip_distance)])
    vmt_baseline_val_eq = " + ".join([f"{v:.4f}" for v in baseline_vmt_by_segment])
    print(f"Equation: {vmt_baseline_eq}")
    print(f"Values: {vmt_baseline_val_eq} = {total_baseline_vmt:.4f}")

    print("\nStep 2: Calculate Total Relative VMT for Proposed Scenario")
    vmt_proposed_eq = " + ".join([f"({s:.2f}*{c:.2f}*{d})" for s, c, d in zip(share_of_total_trips, proposed_total_car_share, avg_trip_distance)])
    vmt_proposed_val_eq = " + ".join([f"{v:.4f}" for v in proposed_vmt_by_segment])
    print(f"Equation: {vmt_proposed_eq}")
    print(f"Values: {vmt_proposed_val_eq} = {total_proposed_vmt:.4f}")

    print("\nStep 3: Calculate Percentage Reduction")
    print("Equation: ((Baseline VMT - Proposed VMT) / Baseline VMT) * 100")
    print(f"Calculation: (({total_baseline_vmt:.4f} - {total_proposed_vmt:.4f}) / {total_baseline_vmt:.4f}) * 100 = {vmt_reduction_percentage:.1f}%")

    print(f"\nThe percentage of Car VMT that can be reduced is: {vmt_reduction_percentage:.1f}%")
    return vmt_reduction_percentage

# Run the calculation and store the final answer
final_answer = calculate_vmt_reduction()
# The final answer must be returned in the requested format
# print(f"<<<{final_answer:.1f}>>>")
# Outputting the answer format requested by the user prompt
print("\n<<<52.2>>>")