import sys
import io

# Redirect stdout to capture print output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

# --- Data from the tables ---
# Common data for both scenarios
share_of_total_trips = [0.48, 0.22, 0.14, 0.16]
avg_trip_distance = [1, 4, 7.5, 20]
distance_categories = ["0-3 miles", "3-5 miles", "5-10 miles", "10+ miles"]

# Baseline data (Table 1)
car_1_occupant_share_baseline = [0.25, 0.38, 0.40, 0.42]
carpool_share_baseline = [0.28, 0.43, 0.45, 0.47]

# Proposed data (Table 2)
car_1_occupant_share_proposed = [0.01, 0.08, 0.14, 0.16]
carpool_share_proposed = [0.02, 0.19, 0.27, 0.33]

# --- Step 1: Calculate Baseline Car VMT ---
total_baseline_vmt = 0
baseline_vmt_parts_str = []
for i in range(len(share_of_total_trips)):
    total_car_share = car_1_occupant_share_baseline[i] + carpool_share_baseline[i]
    vmt_part = share_of_total_trips[i] * total_car_share * avg_trip_distance[i]
    total_baseline_vmt += vmt_part
    # Format the string for the equation part
    part_str = f"({share_of_total_trips[i]} * ({car_1_occupant_share_baseline[i]} + {carpool_share_baseline[i]}) * {avg_trip_distance[i]})"
    baseline_vmt_parts_str.append(part_str)

print("1. Calculate Baseline Car VMT:")
print(f"   Baseline VMT = {' + '.join(baseline_vmt_parts_str)}")
print(f"   Baseline VMT = {total_baseline_vmt:.4f}\n")

# --- Step 2: Calculate Proposed Car VMT ---
total_proposed_vmt = 0
proposed_vmt_parts_str = []
for i in range(len(share_of_total_trips)):
    total_car_share = car_1_occupant_share_proposed[i] + carpool_share_proposed[i]
    vmt_part = share_of_total_trips[i] * total_car_share * avg_trip_distance[i]
    total_proposed_vmt += vmt_part
    # Format the string for the equation part
    part_str = f"({share_of_total_trips[i]} * ({car_1_occupant_share_proposed[i]} + {carpool_share_proposed[i]}) * {avg_trip_distance[i]})"
    proposed_vmt_parts_str.append(part_str)

print("2. Calculate Proposed Car VMT:")
print(f"   Proposed VMT = {' + '.join(proposed_vmt_parts_str)}")
print(f"   Proposed VMT = {total_proposed_vmt:.4f}\n")

# --- Step 3: Calculate Percentage Reduction ---
percentage_reduction = ((total_baseline_vmt - total_proposed_vmt) / total_baseline_vmt) * 100

print("3. Calculate Percentage Reduction:")
print(f"   Percentage Reduction = (({total_baseline_vmt:.4f} - {total_proposed_vmt:.4f}) / {total_baseline_vmt:.4f}) * 100")
print(f"\nFinal Answer: The percentage of Car VMT that can be reduced is {percentage_reduction:.1f}%")

# Restore stdout and print the captured output
sys.stdout = old_stdout
print(captured_output.getvalue())
print(f'<<<{percentage_reduction:.1f}>>>')