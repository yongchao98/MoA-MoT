import sys
import io

# Redirect stdout to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

# --- Data from the tables ---
share_of_total_trips = [0.48, 0.22, 0.14, 0.16]
average_trip_distance = [1, 4, 7.5, 20]
distance_categories = ["0-3 miles", "3-5 miles", "5-10 miles", "10+ miles"]

# Baseline car mode shares (Table 1)
car_1_occupant_baseline = [0.25, 0.38, 0.40, 0.42]
car_pool_baseline = [0.28, 0.43, 0.45, 0.47]

# Proposed car mode shares (Table 2)
car_1_occupant_proposed = [0.01, 0.08, 0.14, 0.16]
car_pool_proposed = [0.02, 0.19, 0.27, 0.33]

# --- Calculation for Baseline VMT ---
vmt_baseline = 0
vmt_baseline_components = []

print("Step 1: Calculate Baseline Car VMT (Vehicle Miles Traveled)")
print("Formula per category: (Share of Trips) * (Total Car Mode Share) * (Avg. Trip Distance)\n")

for i in range(len(distance_categories)):
    total_car_share = car_1_occupant_baseline[i] + car_pool_baseline[i]
    component = share_of_total_trips[i] * total_car_share * average_trip_distance[i]
    vmt_baseline += component
    vmt_baseline_components.append(component)

print(f"Total Baseline VMT = ({share_of_total_trips[0]} * ({car_1_occupant_baseline[0]} + {car_pool_baseline[0]}) * {average_trip_distance[0]}) + "
      f"({share_of_total_trips[1]} * ({car_1_occupant_baseline[1]} + {car_pool_baseline[1]}) * {average_trip_distance[1]}) + "
      f"({share_of_total_trips[2]} * ({car_1_occupant_baseline[2]} + {car_pool_baseline[2]}) * {average_trip_distance[2]}) + "
      f"({share_of_total_trips[3]} * ({car_1_occupant_baseline[3]} + {car_pool_baseline[3]}) * {average_trip_distance[3]})")

print(f"Total Baseline VMT = {vmt_baseline_components[0]:.4f} + {vmt_baseline_components[1]:.4f} + {vmt_baseline_components[2]:.4f} + {vmt_baseline_components[3]:.4f} = {vmt_baseline:.4f}\n")


# --- Calculation for Proposed VMT ---
vmt_proposed = 0
vmt_proposed_components = []

print("Step 2: Calculate Proposed Car VMT")
print("Formula per category: (Share of Trips) * (Total Car Mode Share) * (Avg. Trip Distance)\n")

for i in range(len(distance_categories)):
    total_car_share = car_1_occupant_proposed[i] + car_pool_proposed[i]
    component = share_of_total_trips[i] * total_car_share * average_trip_distance[i]
    vmt_proposed += component
    vmt_proposed_components.append(component)

print(f"Total Proposed VMT = ({share_of_total_trips[0]} * ({car_1_occupant_proposed[0]} + {car_pool_proposed[0]}) * {average_trip_distance[0]}) + "
      f"({share_of_total_trips[1]} * ({car_1_occupant_proposed[1]} + {car_pool_proposed[1]}) * {average_trip_distance[1]}) + "
      f"({share_of_total_trips[2]} * ({car_1_occupant_proposed[2]} + {car_pool_proposed[2]}) * {average_trip_distance[2]}) + "
      f"({share_of_total_trips[3]} * ({car_1_occupant_proposed[3]} + {car_pool_proposed[3]}) * {average_trip_distance[3]})")

print(f"Total Proposed VMT = {vmt_proposed_components[0]:.4f} + {vmt_proposed_components[1]:.4f} + {vmt_proposed_components[2]:.4f} + {vmt_proposed_components[3]:.4f} = {vmt_proposed:.4f}\n")


# --- Calculation for Percentage Reduction ---
reduction = vmt_baseline - vmt_proposed
percentage_reduction = (reduction / vmt_baseline) * 100

print("Step 3: Calculate the Percentage Reduction")
print("Formula: ((Baseline VMT - Proposed VMT) / Baseline VMT) * 100\n")
print(f"Percentage Reduction = (({vmt_baseline:.4f} - {vmt_proposed:.4f}) / {vmt_baseline:.4f}) * 100")
print(f"Final Answer: The percentage of Car VMT that can be reduced is {percentage_reduction:.1f}%")

# Restore stdout
sys.stdout = old_stdout
# Print the captured output
print(captured_output.getvalue())
print(f"<<<{percentage_reduction:.1f}>>>")