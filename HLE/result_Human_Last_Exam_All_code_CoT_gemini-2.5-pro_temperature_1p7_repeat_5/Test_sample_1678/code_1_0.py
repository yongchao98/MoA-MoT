import math

# --- Given Parameters ---
s_hr = 2160.0  # Saturation flow rate (veh/hour)
v_hr = 600.0   # Average flow rate (veh/hour)
R_disp = 56.0  # Displayed red time (s)
Y = 3.0        # Displayed yellow time (s)
AR = 2.0       # All-red time (s)
g = 30.0       # Effective green time (s)
tL = 4.0       # Total lost time (s)
pct_arr_green = 0.40 # Percentage of arrivals during green
pct_arr_red = 0.60   # Percentage of arrivals during red

# --- Plan Execution ---

print("[Step 1] Calculating Cycle Length and Signal Timings")

# Calculate Displayed Green Time (G_disp) using the effective green time formula
# g = G_disp + Y + AR - tL  =>  G_disp = g - Y - AR + tL
G_disp = g - Y - AR + tL
print(f"Based on g = G_disp + Y + AR - tL => {g:.2f} = G_disp + {Y:.2f} + {AR:.2f} - {tL:.2f}")
print(f"Calculated Displayed Green Time (G_disp): {G_disp:.2f} s")

# Calculate Cycle Length (C)
# C = G_disp + Y + R_disp
C = G_disp + Y + R_disp
print(f"Cycle Length (C) = G_disp + Y + R_disp = {G_disp:.2f} + {Y:.2f} + {R_disp:.2f} = {C:.2f} s")

# Calculate Effective Red Time (r)
# r = C - g
r = C - g
print(f"Effective Red Time (r) = C - g = {C:.2f} - {g:.2f} = {r:.2f} s")
print("-" * 40)

print("[Step 2] Calculating Arrival and Service Rates")

# Convert flow rates from veh/hour to veh/second
s = s_hr / 3600.0
v = v_hr / 3600.0
print(f"Saturation Flow Rate (s): {s_hr:.2f} veh/hr = {s:.4f} veh/s")
print(f"Average Arrival Rate (v): {v_hr:.2f} veh/hr = {v:.4f} veh/s")

# Calculate total arrivals per cycle (N_arr)
N_arr = v * C
print(f"Total vehicles per cycle (N_arr) = v * C = {v:.4f} * {C:.2f} = {N_arr:.2f} veh")

# Calculate non-uniform arrival rates (lambda_1 and lambda_2)
lambda_1 = (N_arr * pct_arr_green) / g
lambda_2 = (N_arr * pct_arr_red) / r
print(f"Arrival Rate during Green (λ1) = (N_arr * {pct_arr_green:.2f}) / g = ({N_arr:.2f} * {pct_arr_green:.2f}) / {g:.2f} = {lambda_1:.4f} veh/s")
print(f"Arrival Rate during Red (λ2) = (N_arr * {pct_arr_red:.2f}) / r = ({N_arr:.2f} * {pct_arr_red:.2f}) / {r:.2f} = {lambda_2:.4f} veh/s")
print("-" * 40)


print("[Step 3] Analyzing Queueing Dynamics (D/D/1)")

# Calculate the queue at the end of the effective red period (Q_r)
Q_r = lambda_2 * r
print(f"Queue at end of red (Q_r) = λ2 * r = {lambda_2:.4f} * {r:.2f} = {Q_r:.2f} veh")

# Calculate the time it takes for the queue to clear during green (t_clear)
queue_dissipation_rate = s - lambda_1
t_clear = Q_r / queue_dissipation_rate
print(f"Queue dissipation rate = s - λ1 = {s:.4f} - {lambda_1:.4f} = {queue_dissipation_rate:.4f} veh/s")
print(f"Time to clear queue (t_clear) = Q_r / (s - λ1) = {Q_r:.2f} / {queue_dissipation_rate:.4f} = {t_clear:.2f} s")
print("-" * 40)

print("[Step 4] Calculating Total and Average Delay")

# Calculate the total delay per cycle, which is the sum of two triangle areas in the queue diagram
delay_during_red = 0.5 * r * Q_r
delay_during_green = 0.5 * t_clear * Q_r
total_delay_per_cycle = delay_during_red + delay_during_green

print(f"Delay during red buildup = 0.5 * r * Q_r = 0.5 * {r:.2f} * {Q_r:.2f} = {delay_during_red:.2f} veh-s")
print(f"Delay during green dissipation = 0.5 * t_clear * Q_r = 0.5 * {t_clear:.2f} * {Q_r:.2f} = {delay_during_green:.2f} veh-s")
print(f"Total delay per cycle = {delay_during_red:.2f} + {delay_during_green:.2f} = {total_delay_per_cycle:.2f} veh-s")

# Calculate the average delay per vehicle (d)
avg_delay_per_vehicle = total_delay_per_cycle / N_arr
print(f"Average delay per vehicle (d) = Total Delay / N_arr = {total_delay_per_cycle:.2f} / {N_arr:.2f} = {avg_delay_per_vehicle:.2f} s")
print("-" * 40)

print("\nFinal Answer:")
print(f"The average deterministic delay per vehicle is {avg_delay_per_vehicle:.2f} seconds.")
print(f'<<<{avg_delay_per_vehicle:.2f}>>>')
