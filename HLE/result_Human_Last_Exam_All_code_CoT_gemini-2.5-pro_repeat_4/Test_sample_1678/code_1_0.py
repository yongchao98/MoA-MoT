import math

# Step 1: Define given variables and establish base parameters
s_hr = 2160  # Saturation flow rate (veh/hour)
v_hr = 600   # Approach average flow rate (veh/hour)
g = 30       # Effective green time (s)
R = 56       # Displayed red time (s)
Y = 3        # Displayed yellow time (s)
AR = 2       # All-red time (s)
tL = 4       # Total lost time (s)

# Convert rates to vehicles per second for consistency
s = s_hr / 3600
v = v_hr / 3600

# Calculate the displayed green time (G_disp)
# g = G_disp + Y + AR - tL
G_disp = g - Y - AR + tL

# Calculate the cycle length (C)
# C = G_disp + Y + R
C = G_disp + Y + R

# Calculate the effective red time (r)
# r = C - g
r = C - g

# Step 2: Calculate non-uniform arrival rates (λ1, λ2)
# Total vehicles arriving per cycle (N)
N = v * C

# Number of vehicles arriving during effective green (N_g) and red (N_r)
N_g = 0.40 * N
N_r = 0.60 * N

# Arrival rate during effective green (λ1) and effective red (λ2)
lambda1 = N_g / g
lambda2 = N_r / r

# Step 3: Calculate total delay per cycle (D_cycle)
# The total delay is the area between the cumulative arrival and departure curves.
# It consists of two parts:
# 1. Delay during the red interval: Area1 = 0.5 * lambda2 * r**2
# 2. Delay during the green interval until queue clears: Area2 = 0.5 * (lambda2 * r)**2 / (s - lambda1)
# Note: (lambda2 * r) is the queue length at the end of red, which is N_r.
queue_at_end_of_red = lambda2 * r
delay_during_red = 0.5 * queue_at_end_of_red * r
delay_during_green = 0.5 * (queue_at_end_of_red)**2 / (s - lambda1)
D_cycle = delay_during_red + delay_during_green

# Step 4: Calculate the average deterministic delay per vehicle (d_avg)
d_avg = D_cycle / N

# Print the results and the final equation
print("--- Intermediate Calculations ---")
print(f"Saturation flow rate (s): {s:.4f} veh/s")
print(f"Average arrival rate (v): {v:.4f} veh/s")
print(f"Calculated Displayed Green (G_disp): {G_disp:.2f} s")
print(f"Cycle Length (C): {C:.2f} s")
print(f"Effective Red time (r): {r:.2f} s")
print(f"Arrival rate during green (λ1): {lambda1:.4f} veh/s")
print(f"Arrival rate during red (λ2): {lambda2:.4f} veh/s")
print(f"Total vehicles per cycle (N): {N:.4f} veh")
print(f"Total delay per cycle (D_cycle): {D_cycle:.2f} veh-s")
print("\n--- Final Delay Calculation ---")
print("Average Delay Formula: d = (Total Delay per Cycle) / (Total Vehicles per Cycle)")
print("d = (0.5 * λ2 * r² + 0.5 * (λ2 * r)² / (s - λ1)) / (v * C)")
print(f"d = (0.5 * {lambda2:.4f} * {r:.2f}² + 0.5 * ({lambda2:.4f} * {r:.2f})² / ({s:.4f} - {lambda1:.4f})) / ({v:.4f} * {C:.2f})")
print(f"d = ({delay_during_red:.2f} + {delay_during_green:.2f}) / {N:.2f}")
print(f"d = {D_cycle:.2f} / {N:.2f}")
print(f"\nThe average deterministic delay per vehicle is {d_avg:.2f} seconds.")

<<<23.93>>>