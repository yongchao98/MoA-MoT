import pandas as pd

def calculate_vmt_reduction():
    """
    Calculates the percentage reduction in Car VMT based on the provided tables.
    """
    # Common data for both scenarios
    share_of_total_trips = {
        "0-3 miles": 0.48,
        "3-5 miles": 0.22,
        "5-10 miles": 0.14,
        "10+ miles": 0.16
    }

    average_trip_distance = {
        "0-3 miles": 1,
        "3-5 miles": 4,
        "5-10 miles": 7.5,
        "10+ miles": 20
    }

    # Baseline scenario data (Table 1)
    baseline_car_share = {
        "0-3 miles": 0.25 + 0.28,
        "3-5 miles": 0.38 + 0.43,
        "5-10 miles": 0.40 + 0.45,
        "10+ miles": 0.42 + 0.47
    }

    # Proposed scenario data (Table 2)
    proposed_car_share = {
        "0-3 miles": 0.01 + 0.02,
        "3-5 miles": 0.08 + 0.19,
        "5-10 miles": 0.14 + 0.27,
        "10+ miles": 0.16 + 0.33
    }

    # --- Step 1: Calculate total VMT for the Baseline scenario ---
    total_vmt_baseline = 0
    print("Calculating Baseline Car VMT (relative):")
    for category in share_of_total_trips:
        share = share_of_total_trips[category]
        car_share = baseline_car_share[category]
        distance = average_trip_distance[category]
        vmt_contribution = share * car_share * distance
        total_vmt_baseline += vmt_contribution
        print(f"  - Category '{category}': {share:.2f} (trip share) * {car_share:.2f} (car share) * {distance:.1f} (miles) = {vmt_contribution:.4f}")
    print(f"Total Baseline Car VMT (relative) = {total_vmt_baseline:.4f}\n")

    # --- Step 2: Calculate total VMT for the Proposed scenario ---
    total_vmt_proposed = 0
    print("Calculating Proposed Car VMT (relative):")
    for category in share_of_total_trips:
        share = share_of_total_trips[category]
        car_share = proposed_car_share[category]
        distance = average_trip_distance[category]
        vmt_contribution = share * car_share * distance
        total_vmt_proposed += vmt_contribution
        print(f"  - Category '{category}': {share:.2f} (trip share) * {car_share:.2f} (car share) * {distance:.1f} (miles) = {vmt_contribution:.4f}")
    print(f"Total Proposed Car VMT (relative) = {total_vmt_proposed:.4f}\n")

    # --- Step 3: Calculate the percentage reduction ---
    vmt_reduction = total_vmt_baseline - total_vmt_proposed
    percentage_reduction = (vmt_reduction / total_vmt_baseline) * 100

    print("Calculating Percentage Reduction:")
    print(f"Reduction Amount = Baseline VMT - Proposed VMT = {total_vmt_baseline:.4f} - {total_vmt_proposed:.4f} = {vmt_reduction:.4f}")
    print("\nFinal Equation:")
    print(f"Percentage Reduction = (Reduction Amount / Baseline VMT) * 100")
    print(f"Percentage Reduction = ({vmt_reduction:.4f} / {total_vmt_baseline:.4f}) * 100 = {percentage_reduction:.1f}%")

    return percentage_reduction

# Execute the calculation and print the result
final_answer = calculate_vmt_reduction()
# The final answer is wrapped as requested
# print(f"\n<<< {final_answer:.1f} >>>") # This line is for the final answer format but let's just print the number
print(f"\nBy changing the mode share, the Car VMT can be reduced by {final_answer:.1f}%.")