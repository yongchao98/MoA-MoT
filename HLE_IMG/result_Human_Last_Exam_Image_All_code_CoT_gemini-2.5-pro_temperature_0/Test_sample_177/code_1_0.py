import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Create a string buffer
string_buffer = io.StringIO()
# Redirect stdout to the string buffer
sys.stdout = string_buffer

# --- Data from the tables ---
# These values are the same for both scenarios
share_of_total_trips = [0.48, 0.22, 0.14, 0.16]
avg_trip_distance = [1, 4, 7.5, 20]
distance_labels = ["0-3 miles", "3-5 miles", "5-10 miles", "10+ miles"]

# --- Baseline Data (Table 1) ---
baseline_car_1_occ_share = [0.25, 0.38, 0.40, 0.42]
baseline_car_pool_share = [0.28, 0.43, 0.45, 0.47]

# --- Proposed Data (Table 2) ---
proposed_car_1_occ_share = [0.01, 0.08, 0.14, 0.16]
proposed_car_pool_share = [0.02, 0.19, 0.27, 0.33]

# --- Step 1: Calculate Baseline Car VMT ---
total_vmt_baseline = 0
print("Step 1: Calculate Baseline Car VMT")
print("="*40)
print("Formula: Sum(Share of Total Trips * Total Car Share * Avg Trip Distance)")

# Combine car shares for baseline
baseline_total_car_share = [(s1 + s2) for s1, s2 in zip(baseline_car_1_occ_share, baseline_car_pool_share)]

for i in range(len(share_of_total_trips)):
    vmt_segment = share_of_total_trips[i] * baseline_total_car_share[i] * avg_trip_distance[i]
    total_vmt_baseline += vmt_segment
    print(f"VMT for {distance_labels[i]:<11}: {share_of_total_trips[i]:.2f} * ({baseline_car_1_occ_share[i]:.2f} + {baseline_car_pool_share[i]:.2f}) * {avg_trip_distance[i]:<4} = {vmt_segment:.4f}")

print(f"\nTotal Baseline Car VMT (relative value): {total_vmt_baseline:.4f}")
print("\n")


# --- Step 2: Calculate Proposed Car VMT ---
total_vmt_proposed = 0
print("Step 2: Calculate Proposed Car VMT")
print("="*40)
print("Formula: Sum(Share of Total Trips * Total Car Share * Avg Trip Distance)")

# Combine car shares for proposed
proposed_total_car_share = [(s1 + s2) for s1, s2 in zip(proposed_car_1_occ_share, proposed_car_pool_share)]

for i in range(len(share_of_total_trips)):
    vmt_segment = share_of_total_trips[i] * proposed_total_car_share[i] * avg_trip_distance[i]
    total_vmt_proposed += vmt_segment
    print(f"VMT for {distance_labels[i]:<11}: {share_of_total_trips[i]:.2f} * ({proposed_car_1_occ_share[i]:.2f} + {proposed_car_pool_share[i]:.2f}) * {avg_trip_distance[i]:<4} = {vmt_segment:.4f}")

print(f"\nTotal Proposed Car VMT (relative value): {total_vmt_proposed:.4f}")
print("\n")


# --- Step 3: Calculate Percentage Reduction ---
vmt_reduction = total_vmt_baseline - total_vmt_proposed
percentage_reduction = (vmt_reduction / total_vmt_baseline) * 100

print("Step 3: Calculate Percentage Reduction")
print("="*40)
print(f"VMT Reduction = Baseline VMT - Proposed VMT")
print(f"              = {total_vmt_baseline:.4f} - {total_vmt_proposed:.4f} = {vmt_reduction:.4f}")
print(f"\nPercentage Reduction = (VMT Reduction / Baseline VMT) * 100")
print(f"                     = ({vmt_reduction:.4f} / {total_vmt_baseline:.4f}) * 100")
print(f"                     = {percentage_reduction:.1f}%")

print(f"\nBy changing the mode share, the Car VMT can be reduced by {percentage_reduction:.1f}%.")

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = string_buffer.getvalue()
# Print the output to the real console
print(output)