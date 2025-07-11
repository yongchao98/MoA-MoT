import math

# Step 1: Define given parameters
s_h = 2160.0  # Saturation flow rate (veh/hour)
R = 56.0      # Displayed red time (seconds)
Y = 3.0       # Displayed yellow time (seconds)
AR = 2.0      # All-red time (seconds)
g = 30.0      # Effective green time (seconds)
t_L = 4.0     # Total lost time (seconds)
v_h = 600.0   # Approach average flow rate (veh/hour)
pct_green_arrival = 0.40 # Percentage of traffic arriving during green
pct_red_arrival = 0.60   # Percentage of traffic arriving during red

# Step 2: Determine Cycle Timing
# The formula for effective green time is g = G_d + Y + AR - t_L. We can find G_d from it.
G_d = g - Y - AR + t_L
# The cycle length C is the sum of all displayed intervals for the approach.
C = G_d + Y + R
# The effective red time r is the portion of the cycle that is not effectively green.
r = C - g

# Step 3: Convert flow rates to vehicles per second
s = s_h / 3600.0
v = v_h / 3600.0

# Step 4: Calculate piecewise arrival rates (λ1 and λ2)
# Total vehicles arriving in one cycle (N)
N = v * C
# Arrivals during the green interval (N_g) and red interval (N_r)
N_g = N * pct_green_arrival
N_r = N * pct_red_arrival
# Rate of arrival during green (λ1) and red (λ2)
lambda1 = N_g / g
lambda2 = N_r / r

# Step 5: Calculate total delay per cycle (D_total)
# The total delay is the area under the queue length curve.
# It can be calculated as the sum of delay during red and delay until the queue dissipates during green.
# Delay during red = integral of queue length from 0 to r = 0.5 * λ2 * r^2
delay_during_red_phase = 0.5 * lambda2 * r**2
# Queue at the end of red = λ2 * r
queue_at_end_of_red = lambda2 * r
# Delay during green phase until queue clears = Area of a triangle with height (queue_at_end_of_red)
# and base (time to clear queue) = (queue_at_end_of_red)^2 / (2 * (s - λ1))
delay_during_green_phase = (queue_at_end_of_red**2) / (2 * (s - lambda1))
D_total = delay_during_red_phase + delay_during_green_phase

# Step 6: Calculate average delay per vehicle (d)
# Average delay is the total delay divided by the total number of vehicles in a cycle.
d = D_total / N

# Print the results and the final equation with numbers
print("--- Calculation of Cycle Timing ---")
print(f"Displayed Green Time (G_d) = {g:.2f}s - {Y:.2f}s - {AR:.2f}s + {t_L:.2f}s = {G_d:.2f} s")
print(f"Cycle Length (C) = {G_d:.2f}s + {Y:.2f}s + {R:.2f}s = {C:.2f} s")
print(f"Effective Red Time (r) = {C:.2f}s - {g:.2f}s = {r:.2f} s")

print("\n--- Calculation of Flow Rates (per second) ---")
print(f"Saturation Flow Rate (s) = {s_h:.2f} veh/h / 3600 = {s:.4f} veh/s")
print(f"Average Flow Rate (v) = {v_h:.2f} veh/h / 3600 = {v:.4f} veh/s")
print(f"Arrival Rate during Green (λ1) = ({v:.4f} * {C:.2f} * {pct_green_arrival:.2f}) / {g:.2f} = {lambda1:.4f} veh/s")
print(f"Arrival Rate during Red (λ2) = ({v:.4f} * {C:.2f} * {pct_red_arrival:.2f}) / {r:.2f} = {lambda2:.4f} veh/s")

print("\n--- Final Average Delay Calculation ---")
print("Average Delay (d) = Total Delay per Cycle (D_total) / Total Vehicles per Cycle (N)")
print(f"D_total = (0.5 * λ2 * r^2) + ((λ2 * r)^2 / (2 * (s - λ1)))")
print(f"D_total = (0.5 * {lambda2:.4f} * {r:.2f}^2) + (({lambda2:.4f} * {r:.2f})^2 / (2 * ({s:.4f} - {lambda1:.4f})))")
print(f"D_total = {delay_during_red_phase:.2f} veh-s + {delay_during_green_phase:.2f} veh-s = {D_total:.2f} veh-s")
print(f"N = {v:.4f} * {C:.2f} = {N:.2f} vehicles")
print(f"d = {D_total:.2f} / {N:.2f}")

print(f"\nThe average deterministic delay per vehicle is {d:.2f} seconds.")
<<<23.93>>>