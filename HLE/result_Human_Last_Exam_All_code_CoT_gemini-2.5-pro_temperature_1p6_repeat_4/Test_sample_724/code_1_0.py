import math

# --- Parameters ---
# n: The number of 100-liter containers of water at the start.
n = 3
# m: The distance to the destination in kilometers.
m = 50
# ------------------

# Total initial water
initial_water = n * 100.0

# Find k, the number of full depot-establishing legs
k = 0
traveled_dist_sum = 0.0
for i in range(1, n):
    # Number of containers at the start of this leg is n - i + 1
    # Number of trips is 2 * (n - i + 1) - 1
    num_trips = 2 * n - 2 * i + 1
    leg_dist = 100.0 / num_trips
    
    if m >= traveled_dist_sum + leg_dist:
        traveled_dist_sum += leg_dist
        k = i
    else:
        # Cannot complete this leg, so break
        break

# --- Calculate the final result ---
sum_of_leg_distances = 0.0
for i in range(1, k + 1):
    num_trips = 2 * n - 2 * i + 1
    sum_of_leg_distances += 100.0 / num_trips

remaining_distance = m - sum_of_leg_distances
containers_on_last_leg = n - k
trips_on_last_leg = 2 * containers_on_last_leg - 1

water_consumed_on_full_legs = k * 100.0
water_consumed_on_last_leg = remaining_distance * trips_on_last_leg

total_water_consumed = water_consumed_on_full_legs + water_consumed_on_last_leg
water_left = initial_water - total_water_consumed

# --- Construct the output equation string ---
sum_terms = []
for i in range(1, k + 1):
    denominator = 2 * n - 2 * i + 1
    sum_terms.append(f"100/{denominator}")

sum_str = " + ".join(sum_terms) if sum_terms else "0"

# Expression for the amount of water remaining at the start of the final leg
water_at_start_of_last_leg_expr = f"({n}-{k})*100"

# Expression for the water consumed on the final leg
water_consumed_on_last_leg_expr = f"({m} - ({sum_str})) * {trips_on_last_leg}"

# Final equation
final_equation = f"{water_at_start_of_last_leg_expr} - {water_consumed_on_last_leg_expr}"

print(f"For n={n} and m={m} km:")
print(f"The maximum amount of water left is: {water_left:.4f} liters")
print("\nThe calculation is based on the formula: Water Left = Initial Water - Consumed Water")
print("Expressed in terms of depot stages:")
print(f"Water Left = {final_equation}")
