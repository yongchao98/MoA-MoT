# Data from Table 1 (Baseline) and Table 2 (Proposed)
share_of_total_trips = [0.48, 0.22, 0.14, 0.16]
avg_trip_distances = [1, 4, 7.5, 20]
trip_distance_labels = ["0-3 miles", "3-5 miles", "5-10 miles", "10+ miles"]

# Baseline car mode shares
baseline_car1_share = [0.25, 0.38, 0.40, 0.42]
baseline_car2_share = [0.28, 0.43, 0.45, 0.47]

# Proposed car mode shares
proposed_car1_share = [0.01, 0.08, 0.14, 0.16]
proposed_car2_share = [0.02, 0.19, 0.27, 0.33]

# --- 1. Calculate Baseline Car VMT ---
baseline_vmt = 0
baseline_vmt_calcs = []

print("Calculating Baseline Car VMT...")
for i in range(len(share_of_total_trips)):
    total_car_share = baseline_car1_share[i] + baseline_car2_share[i]
    vmt_component = share_of_total_trips[i] * total_car_share * avg_trip_distances[i]
    baseline_vmt += vmt_component
    print(f"  VMT for {trip_distance_labels[i]:>10}: {share_of_total_trips[i]} (trip share) * ({baseline_car1_share[i]} + {baseline_car2_share[i]}) (car mode share) * {avg_trip_distances[i]} (miles) = {vmt_component:.4f}")
    baseline_vmt_calcs.append(f"{vmt_component:.4f}")

print(f"\nTotal Baseline Car VMT = {' + '.join(baseline_vmt_calcs)} = {baseline_vmt:.4f}\n")


# --- 2. Calculate Proposed Car VMT ---
proposed_vmt = 0
proposed_vmt_calcs = []

print("Calculating Proposed Car VMT...")
for i in range(len(share_of_total_trips)):
    total_car_share = proposed_car1_share[i] + proposed_car2_share[i]
    vmt_component = share_of_total_trips[i] * total_car_share * avg_trip_distances[i]
    proposed_vmt += vmt_component
    print(f"  VMT for {trip_distance_labels[i]:>10}: {share_of_total_trips[i]} (trip share) * ({proposed_car1_share[i]} + {proposed_car2_share[i]}) (car mode share) * {avg_trip_distances[i]} (miles) = {vmt_component:.4f}")
    proposed_vmt_calcs.append(f"{vmt_component:.4f}")

print(f"\nTotal Proposed Car VMT = {' + '.join(proposed_vmt_calcs)} = {proposed_vmt:.4f}\n")


# --- 3. Calculate Percentage Reduction ---
vmt_reduction = baseline_vmt - proposed_vmt
percentage_reduction = (vmt_reduction / baseline_vmt) * 100

print("Calculating Percentage Reduction...")
print(f"Percentage Reduction = (({baseline_vmt:.4f} - {proposed_vmt:.4f}) / {baseline_vmt:.4f}) * 100")
print(f"Result: {percentage_reduction:.1f}%")

# Final Answer
print(f"\nBy changing the mode share, Car VMT can be reduced by {percentage_reduction:.1f}%.")
print(f"<<<{percentage_reduction:.1f}>>>")