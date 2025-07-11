import math

# Given parameters
s_hr = 2160  # Saturation flow rate (veh/hour)
R = 56       # Displayed red time (s)
Y = 3        # Displayed yellow time (s)
AR = 2       # All-red time (s)
g = 30       # Effective green time (s)
t_L = 4      # Total lost time (s)
v_hr = 600   # Approach average flow rate (veh/hour)
percent_green_arrival = 0.40 # 40% of traffic arrives during green
percent_red_arrival = 0.60   # 60% of traffic arrives during red

# Step 1: Establish Cycle Parameters
# Effective red time (r) is the sum of displayed red and all-red time.
r = R + AR
# Cycle length (C) is the sum of effective green and effective red time.
C = g + r

print(f"Step 1: Cycle Parameters")
print(f"Effective red time (r) = {R}s + {AR}s = {r}s")
print(f"Cycle length (C) = {g}s + {r}s = {C}s\n")

# Step 2: Convert Units
# Convert flow rates from veh/hour to veh/s
s = s_hr / 3600
v = v_hr / 3600

print(f"Step 2: Flow Rates in veh/s")
print(f"Saturation flow rate (s) = {s_hr}/3600 = {s:.4f} veh/s")
print(f"Average arrival rate (v) = {v_hr}/3600 = {v:.4f} veh/s\n")

# Step 3: Calculate Arrival Numbers and Rates
# Total vehicles per cycle (N)
N = v * C
# Number of vehicles arriving during green (N_g) and red (N_r)
N_g = percent_green_arrival * N
N_r = percent_red_arrival * N
# Arrival rate during green (lambda1) and red (lambda2)
lambda1 = N_g / g
lambda2 = N_r / r

print(f"Step 3: Arrival Calculations")
print(f"Total vehicles per cycle (N) = {v:.4f} * {C} = {N:.4f} veh")
print(f"Vehicles arriving during red (N_r) = {percent_red_arrival} * {N:.4f} = {N_r:.4f} veh")
print(f"Vehicles arriving during green (N_g) = {percent_green_arrival} * {N:.4f} = {N_g:.4f} veh")
print(f"Arrival rate during red (λ2) = {N_r:.4f} / {r} = {lambda2:.4f} veh/s")
print(f"Arrival rate during green (λ1) = {N_g:.4f} / {g} = {lambda1:.4f} veh/s\n")

# Step 4: Calculate Total Delay per Cycle (D_total)
# Delay accumulated during the red interval
delay_red_part = 0.5 * lambda2 * r**2
# Delay accumulated during the green interval until queue clears
# Note: (lambda2 * r) is the queue size at the start of green, which is N_r
queue_at_green_start = lambda2 * r
delay_green_part = (0.5 * (queue_at_green_start)**2) / (s - lambda1)
# Total delay per cycle
D_total = delay_red_part + delay_green_part

print(f"Step 4: Total Delay Calculation")
print("Total Delay (D_total) = Delay during red + Delay during green queue")
print(f"D_total = 0.5 * {lambda2:.4f} * {r}**2 + (0.5 * ({lambda2:.4f} * {r})**2) / ({s:.4f} - {lambda1:.4f})")
print(f"D_total = {delay_red_part:.2f} + {delay_green_part:.2f}")
print(f"D_total = {D_total:.2f} veh-s\n")

# Step 5: Calculate Average Delay per Vehicle (d)
# Average delay per vehicle is the total delay divided by the number of vehicles
d = D_total / N

print(f"Step 5: Average Delay per Vehicle")
print(f"Average Delay (d) = Total Delay (D_total) / Total Arrivals (N)")
print(f"d = {D_total:.2f} / {N:.2f}")
print(f"d = {d:.2f} seconds\n")

print("Final Answer:")
print(f"The average deterministic delay per vehicle is {d:.2f} seconds.")
<<<23.93>>>