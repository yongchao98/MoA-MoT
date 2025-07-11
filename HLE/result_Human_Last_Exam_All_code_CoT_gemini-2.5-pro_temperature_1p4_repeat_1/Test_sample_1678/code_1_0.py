import math

# Step 1: Define Inputs
s_hr = 2160  # Saturation flow rate (veh/hour)
R_disp = 56  # Displayed red time (seconds)
Y = 3        # Displayed yellow time (seconds)
AR = 2       # All-red time (seconds)
g = 30       # Effective green time (seconds)
tL = 4       # Total lost time (seconds)
v_hr = 600   # Approach average flow rate (veh/hour)
pct_green_arr = 0.40 # Percentage of traffic arriving during effective green
pct_red_arr = 0.60   # Percentage of traffic arriving during effective red

print("--- Given Parameters ---")
print(f"Saturation flow rate (s): {s_hr} veh/hour")
print(f"Average flow rate (v): {v_hr} veh/hour")
print(f"Effective green time (g): {g} s")
print(f"Displayed red time (R_disp): {R_disp} s")
print(f"Displayed yellow time (Y): {Y} s")
print(f"All-red time (AR): {AR} s")
print(f"Total lost time (tL): {tL} s\n")

# Step 2: Calculate Cycle Parameters
# Total lost time tL = startup lost time (l1) + clearance lost time (l2)
# Assume clearance lost time is the all-red time.
l2 = AR
l1 = tL - l2
# Effective green g = Displayed Green (G) + Yellow (Y) - startup lost time (l1)
# So, G can be calculated.
G = g - Y + l1
# Cycle Length C = Displayed Green (G) + Displayed Yellow (Y) + Displayed Red (R_disp)
C = G + Y + R_disp
# Effective red time r = Cycle Length (C) - Effective Green (g)
r = C - g

print("--- Step 1: Calculating Cycle Parameters ---")
print(f"Startup lost time (l1 = tL - AR): {tL} - {AR} = {l1} s")
print(f"Displayed green time (G = g - Y + l1): {g} - {Y} + {l1} = {G} s")
print(f"Cycle length (C = G + Y + R_disp): {G} + {Y} + {R_disp} = {C} s")
print(f"Effective red time (r = C - g): {C} - {g} = {r} s\n")


# Step 3: Calculate Traffic Flow Parameters
# Convert hourly rates to vehicles per second
s = s_hr / 3600.0
v = v_hr / 3600.0

# Total vehicles arriving per cycle
N_total = v * C

# Number of vehicles arriving during green and red intervals
N_g = N_total * pct_green_arr
N_r = N_total * pct_red_arr

# Arrival rate during effective green (lambda1) and effective red (lambda2)
lambda1 = N_g / g
lambda2 = N_r / r

print("--- Step 2: Calculating Traffic Flow Parameters ---")
print(f"Saturation flow rate (s): {s_hr}/3600 = {s:.4f} veh/s")
print(f"Average flow rate (v): {v_hr}/3600 = {v:.4f} veh/s")
print(f"Total vehicles per cycle (N = v * C): {v:.4f} * {C} = {N_total:.4f} veh")
print(f"Vehicles arriving in red (N_r): {N_total:.4f} * {pct_red_arr} = {N_r:.4f} veh")
print(f"Arrival rate in red (λ2 = N_r / r): {N_r:.4f} / {r} = {lambda2:.4f} veh/s")
print(f"Vehicles arriving in green (N_g): {N_total:.4f} * {pct_green_arr} = {N_g:.4f} veh")
print(f"Arrival rate in green (λ1 = N_g / g): {N_g:.4f} / {g} = {lambda1:.4f} veh/s\n")

# Step 4: Analyze the Queue
# Time to clear the queue that built up during red
t_clear = N_r / (s - lambda1)

print("--- Step 3: Calculating Queue Clearing Time ---")
print(f"Queue at start of green (Q_max = N_r): {N_r:.4f} veh")
print(f"Time to clear queue (t_clear = N_r / (s - λ1)): {N_r:.4f} / ({s:.4f} - {lambda1:.4f}) = {t_clear:.4f} s\n")

# Step 5: Calculate Total and Average Delay
# Total delay per cycle is the area under the queueing diagram.
# It can be calculated as the sum of delay during red and delay during queue clearance in green.
# Delay during red = 0.5 * r * N_r
delay_red_period = 0.5 * r * N_r
# Delay during green period (while queue clears) = 0.5 * t_clear * N_r
delay_green_period = 0.5 * t_clear * N_r
# Total delay per cycle
D_total = delay_red_period + delay_green_period

# Average delay per vehicle
d_avg = D_total / N_total

print("--- Step 4: Calculating Average Delay per Vehicle ---")
print(f"Total delay per cycle (D_total = 0.5*N_r*r + 0.5*N_r*t_clear): 0.5*{N_r:.2f}*{r} + 0.5*{N_r:.2f}*{t_clear:.2f} = {D_total:.2f} veh-s")
print(f"Average delay per vehicle (d = D_total / N_total): {D_total:.2f} / {N_total:.2f} = {d_avg:.2f} s\n")

# Step 6: Final Answer
print("--- Final Answer ---")
print(f"The average deterministic delay per vehicle is {d_avg:.2f} seconds.")
