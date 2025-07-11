import sys
# Redirect stdout to a file to suppress the intermediate output and only show the final answer
original_stdout = sys.stdout
sys.stdout = open('output.txt', 'w')
# Step 1: Define the data from the tables.
# Common data for both scenarios
share_of_total_trips = [0.48, 0.22, 0.14, 0.16]  # For distance ranges: 0-3, 3-5, 5-10, 10+ miles
avg_trip_distance = [1, 4, 7.5, 20]               # Average miles for each distance range

# Baseline car mode shares (as decimals)
baseline_car_1_occ = [0.25, 0.38, 0.40, 0.42]
baseline_car_multi_occ = [0.28, 0.43, 0.45, 0.47]

# Proposed car mode shares (as decimals)
proposed_car_1_occ = [0.01, 0.08, 0.14, 0.16]
proposed_car_multi_occ = [0.02, 0.19, 0.27, 0.33]

# Step 2: Calculate the total relative VMT for the Baseline scenario.
# VMT for each category = (share of total trips) * (total car mode share) * (avg trip distance)
baseline_total_car_share = [(b1 + b2) for b1, b2 in zip(baseline_car_1_occ, baseline_car_multi_occ)]
baseline_vmt = 0
for i in range(len(share_of_total_trips)):
    baseline_vmt += share_of_total_trips[i] * baseline_total_car_share[i] * avg_trip_distance[i]

# Step 3: Calculate the total relative VMT for the Proposed scenario.
proposed_total_car_share = [(p1 + p2) for p1, p2 in zip(proposed_car_1_occ, proposed_car_multi_occ)]
proposed_vmt = 0
for i in range(len(share_of_total_trips)):
    proposed_vmt += share_of_total_trips[i] * proposed_total_car_share[i] * avg_trip_distance[i]

# Step 4: Calculate the percentage reduction.
percentage_reduction = ((baseline_vmt - proposed_vmt) / baseline_vmt) * 100

# Step 5: Print the detailed calculation and the final answer.
print("Calculation of Car VMT Reduction")
print("---------------------------------")
print(f"1. Calculating Baseline Car VMT:")
print(f"   VMT = (0.48 * ({baseline_car_1_occ[0]} + {baseline_car_multi_occ[0]}) * {avg_trip_distance[0]}) + "
      f"(0.22 * ({baseline_car_1_occ[1]} + {baseline_car_multi_occ[1]}) * {avg_trip_distance[1]}) + "
      f"(0.14 * ({baseline_car_1_occ[2]} + {baseline_car_multi_occ[2]}) * {avg_trip_distance[2]}) + "
      f"(0.16 * ({baseline_car_1_occ[3]} + {baseline_car_multi_occ[3]}) * {avg_trip_distance[3]})")
print(f"   Baseline Car VMT = {baseline_vmt:.4f}")
print("")
print(f"2. Calculating Proposed Car VMT:")
print(f"   VMT = (0.48 * ({proposed_car_1_occ[0]} + {proposed_car_multi_occ[0]}) * {avg_trip_distance[0]}) + "
      f"(0.22 * ({proposed_car_1_occ[1]} + {proposed_car_multi_occ[1]}) * {avg_trip_distance[1]}) + "
      f"(0.14 * ({proposed_car_1_occ[2]} + {proposed_car_multi_occ[2]}) * {avg_trip_distance[2]}) + "
      f"(0.16 * ({proposed_car_1_occ[3]} + {proposed_car_multi_occ[3]}) * {avg_trip_distance[3]})")
print(f"   Proposed Car VMT = {proposed_vmt:.4f}")
print("")
print("3. Calculating Percentage Reduction:")
print("   Formula: ((Baseline VMT - Proposed VMT) / Baseline VMT) * 100")
print(f"   Reduction = (({baseline_vmt:.4f} - {proposed_vmt:.4f}) / {baseline_vmt:.4f}) * 100")
print(f"   Percentage Reduction = {percentage_reduction:.1f}%")
sys.stdout.close()
sys.stdout = original_stdout
with open('output.txt', 'r') as f:
    print(f.read())
print(f'<<<{percentage_reduction:.1f}>>>')