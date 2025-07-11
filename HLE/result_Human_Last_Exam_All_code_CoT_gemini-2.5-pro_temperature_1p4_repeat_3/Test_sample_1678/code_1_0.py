import math

# Step 1: Define given parameters
s_hr = 2160  # Saturation flow rate (veh/hour)
R_disp = 56  # Displayed red time (seconds)
Y_disp = 3   # Displayed yellow time (seconds)
AR = 2       # All-red time (seconds)
g = 30       # Effective green time (seconds)
t_L = 4      # Total lost time (seconds)
v_hr = 600   # Average approach flow rate (veh/hour)
perc_green_arrival = 0.40 # Percentage of traffic arriving during green
perc_red_arrival = 0.60   # Percentage of traffic arriving during red

# Step 2: Calculate cycle parameters
# Formula for effective green time: g = G_disp + Y_disp + AR - t_L
# We can find the displayed green time (G_disp) from this.
G_disp = g - Y_disp - AR + t_L
# Formula for cycle length: C = G_disp + Y_disp + R_disp
C = G_disp + Y_disp + R_disp
# Effective red time r = C - g
r = C - g

# Step 3: Convert rates to vehicles per second
s = s_hr / 3600.0
v = v_hr / 3600.0

# Step 4: Calculate arrival rates for green and red intervals
# Total number of vehicles arriving per cycle
N_total = v * C
# Number of vehicles arriving during effective green and red
N_green = N_total * perc_green_arrival
N_red = N_total * perc_red_arrival
# Arrival rate during effective green (lambda1)
lambda1 = N_green / g
# Arrival rate during effective red (lambda2)
lambda2 = N_red / r

# Step 5: Calculate time to clear the queue (t_c) during the green interval
# The queue length at the end of the red interval is Q_r = lambda2 * r
# The queue dissipates at a rate of (s - lambda1)
# t_c = Q_r / (s - lambda1)
t_c = (lambda2 * r) / (s - lambda1)

# Step 6: Calculate total delay per cycle (D_total)
# D_total is the area under the queue vs. time graph.
# It can be calculated as the sum of two areas:
# 1. Area of triangle during red time: 0.5 * r * (lambda2 * r)
# 2. Area of triangle during queue clearance in green time: 0.5 * t_c * (lambda2 * r)
# D_total = 0.5 * (lambda2 * r) * (r + t_c)
D_total = 0.5 * lambda2 * r * (r + t_c)

# Step 7: Calculate average delay per vehicle (d)
d = D_total / N_total

# Step 8: Print the results and the final equation
print("--- Intermediate Calculations ---")
print(f"Displayed Green Time (G_disp): {G_disp:.2f} s")
print(f"Cycle Length (C): {C:.2f} s")
print(f"Effective Red Time (r): {r:.2f} s")
print(f"Saturation Flow Rate (s): {s:.4f} veh/s")
print(f"Average Arrival Rate (v): {v:.4f} veh/s")
print(f"Arrival Rate during Green (λ1): {lambda1:.4f} veh/s")
print(f"Arrival Rate during Red (λ2): {lambda2:.4f} veh/s")
print(f"Time to Clear Queue (t_c): {t_c:.2f} s")
print(f"Total Vehicles per Cycle (N_total): {N_total:.2f} veh")
print(f"Total Delay per Cycle (D_total): {D_total:.2f} veh-s")
print("\n--- Final Calculation for Average Delay per Vehicle ---")
print("Average Delay (d) = Total Delay per Cycle / Total Vehicles per Cycle")
print(f"d = {D_total:.2f} / {N_total:.2f}")
print(f"\nThe average deterministic delay per vehicle is: {d:.2f} seconds.")

# Final answer in the required format
final_answer = round(d, 2)
# <<<23.93>>>