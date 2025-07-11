import numpy as np

def calculate_vmt(car_1_occupant_share, car_multi_occupant_share, share_of_total_trips, average_trip_distance):
    """
    Calculates the total relative Vehicle Miles Traveled (VMT).
    VMT is calculated as the sum of VMT for each trip distance category.
    VMT per category = (share of total trips) * (total car mode share) * (average trip distance)
    """
    # Sum the shares of the two car modes to get the total car mode share for each distance category
    total_car_share = np.add(car_1_occupant_share, car_multi_occupant_share)
    
    # Calculate the VMT for each distance category
    vmt_by_category = np.multiply(np.multiply(share_of_total_trips, total_car_share), average_trip_distance)
    
    # Sum the VMT across all categories to get the total VMT
    total_vmt = np.sum(vmt_by_category)
    
    return total_vmt

# --- Data from the tables ---

# Data that is common for both scenarios
share_of_total_trips = np.array([0.48, 0.22, 0.14, 0.16])
average_trip_distance = np.array([1, 4, 7.5, 20]) # in miles

# Baseline scenario data (Table 1)
baseline_car_1_occupant_share = np.array([0.25, 0.38, 0.40, 0.42])
baseline_car_multi_occupant_share = np.array([0.28, 0.43, 0.45, 0.47])

# Proposed scenario data (Table 2)
proposed_car_1_occupant_share = np.array([0.01, 0.08, 0.14, 0.16])
proposed_car_multi_occupant_share = np.array([0.02, 0.19, 0.27, 0.33])


# --- Calculations ---

# 1. Calculate Baseline Car VMT
baseline_vmt = calculate_vmt(baseline_car_1_occupant_share, baseline_car_multi_occupant_share, share_of_total_trips, average_trip_distance)

# 2. Calculate Proposed Car VMT
proposed_vmt = calculate_vmt(proposed_car_1_occupant_share, proposed_car_multi_occupant_share, share_of_total_trips, average_trip_distance)

# 3. Calculate the percentage reduction
vmt_reduction = baseline_vmt - proposed_vmt
percentage_reduction = (vmt_reduction / baseline_vmt) * 100

# --- Output the results ---

print("--- Calculating Car VMT Reduction ---")
print(f"Relative Baseline Car VMT = {baseline_vmt:.4f}")
print(f"Relative Proposed Car VMT = {proposed_vmt:.4f}\n")

print("--- Final Equation ---")
print(f"Percentage Reduction = ( (Baseline VMT - Proposed VMT) / Baseline VMT ) * 100")
print(f"Percentage Reduction = ( ({baseline_vmt:.4f} - {proposed_vmt:.4f}) / {baseline_vmt:.4f} ) * 100")
print(f"Percentage Reduction = ( {vmt_reduction:.4f} / {baseline_vmt:.4f} ) * 100")
print(f"Result = {percentage_reduction:.1f}%")

# The final numerical answer
print("\nThe percentage of Car VMT that can be reduced is:")
print(f"{percentage_reduction:.1f}%")
<<<52.2>>>