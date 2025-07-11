import sys
import io

# Redirect stdout to capture print output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

# --- Data from the tables ---

# Share of total trips and average distance for each category (same for both tables)
# Categories: 0-3 miles, 3-5 miles, 5-10 miles, 10+ miles
share_of_total_trips = [0.48, 0.22, 0.14, 0.16]
avg_trip_distance = [1, 4, 7.5, 20]

# --- Baseline Data (Table 1) ---
# Mode share for 'car, 1 occupant'
baseline_car1_share = [0.25, 0.38, 0.40, 0.42]
# Mode share for 'car, with more than 2 occupants'
baseline_car2_share = [0.28, 0.43, 0.45, 0.47]

# --- Proposed Data (Table 2) ---
# Mode share for 'car, 1 occupant'
proposed_car1_share = [0.01, 0.08, 0.14, 0.16]
# Mode share for 'car, with more than 2 occupants'
proposed_car2_share = [0.02, 0.19, 0.27, 0.33]

# --- Calculation ---

def calculate_vmt(car1_shares, car2_shares):
    """Calculates the total relative VMT based on car mode shares."""
    total_vmt = 0
    for i in range(len(share_of_total_trips)):
        # Total car mode share for the distance category
        total_car_share = car1_shares[i] + car2_shares[i]
        # Weighted VMT for the category
        vmt_category = share_of_total_trips[i] * total_car_share * avg_trip_distance[i]
        total_vmt += vmt_category
    return total_vmt

# Calculate VMT for both scenarios
baseline_vmt = calculate_vmt(baseline_car1_share, baseline_car2_share)
proposed_vmt = calculate_vmt(proposed_car1_share, proposed_car2_share)

# Calculate the percentage reduction
percentage_reduction = ((baseline_vmt - proposed_vmt) / baseline_vmt) * 100

# --- Output ---

print(f"Baseline Car VMT: {baseline_vmt:.4f}")
print(f"Proposed Car VMT: {proposed_vmt:.4f}")
print("\nTo find the percentage reduction in Car VMT, we use the formula:")
print("Percentage Reduction = ((Baseline VMT - Proposed VMT) / Baseline VMT) * 100")
print("\nPlugging in the calculated values:")
print(f"Percentage Reduction = (({baseline_vmt:.4f} - {proposed_vmt:.4f}) / {baseline_vmt:.4f}) * 100")
print(f"Percentage Reduction = {baseline_vmt - proposed_vmt:.4f} / {baseline_vmt:.4f} * 100")
print(f"The percentage of Car VMT that can be reduced is: {percentage_reduction:.1f}%")

# Restore stdout
sys.stdout = old_stdout
# Get the captured output
output = captured_output.getvalue()
print(output)

# Final answer in the required format
print(f"<<<{percentage_reduction:.1f}>>>")