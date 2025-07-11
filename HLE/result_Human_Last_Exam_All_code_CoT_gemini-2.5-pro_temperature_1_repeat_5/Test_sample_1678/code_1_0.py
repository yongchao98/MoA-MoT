import math

# Step 1: Define given parameters
s_h = 2160  # Saturation flow rate (veh/hour)
v_h = 600   # Average approach flow rate (veh/hour)
R_disp = 56 # Displayed red time (seconds)
Y = 3       # Displayed yellow time (seconds)
AR = 2      # All-red time (seconds)
g = 30      # Effective green time (seconds)
# t_L = 4     # Total lost time (seconds) - Not directly needed for this calculation method

# Convert flow rates to veh/second
s = s_h / 3600
v = v_h / 3600

print("--- Step 1: Initial Parameters ---")
print(f"Saturation flow rate (s): {s_h} veh/h or {s:.4f} veh/s")
print(f"Average flow rate (v): {v_h} veh/h or {v:.4f} veh/s")
print(f"Effective green time (g): {g} s")

# Step 2: Calculate cycle timings
# The effective red time (r) is the period when the queue forms.
# It is the sum of displayed red, yellow, and all-red.
r = R_disp + Y + AR
# The cycle length (C) is the sum of effective green and effective red times.
C = g + r

print("\n--- Step 2: Cycle Timings ---")
print(f"Effective red time (r) = {R_disp} + {Y} + {AR} = {r} s")
print(f"Cycle length (C) = {g} + {r} = {C} s")

# Step 3: Calculate arrival rates
# Total number of vehicles arriving per cycle
N = v * C
# 60% of traffic arrives during the effective red interval (r)
vehicles_during_red = 0.60 * N
# 40% of traffic arrives during the effective green interval (g)
vehicles_during_green = 0.40 * N

# Arrival rate during red (lambda_2)
lambda_2 = vehicles_during_red / r
# Arrival rate during green (lambda_1)
lambda_1 = vehicles_during_green / g

print("\n--- Step 3: Arrival Rates ---")
print(f"Total vehicles per cycle (N): {N:.4f}")
print(f"Arrival rate during red (λ2): {lambda_2:.4f} veh/s")
print(f"Arrival rate during green (λ1): {lambda_1:.4f} veh/s")

# Step 4: Analyze the queue
# Queue at the start of green (Q₀) is the number of vehicles that arrived during red.
Q0 = lambda_2 * r

# Time to clear the queue (t_c)
# The queue dissipates at a rate of (s - lambda_1)
t_c = Q0 / (s - lambda_1)

print("\n--- Step 4: Queue Analysis ---")
print(f"Queue at start of green (Q₀): {Q0:.4f} vehicles")
print(f"Time to clear queue (t_c): {t_c:.4f} s")
if t_c > g:
    print("Warning: The queue does not clear within the effective green time.")

# Step 5: Calculate total delay per cycle
# Total delay (D_total) is the area of the queueing diagram.
# This area can be calculated as the area of a polygon: 0.5 * Q₀ * (r + t_c)
D_total = 0.5 * Q0 * (r + t_c)

print("\n--- Step 5: Total Delay Calculation ---")
print(f"Total delay per cycle (D_total) = 0.5 * {Q0:.2f} * ({r} + {t_c:.2f}) = {D_total:.2f} veh-seconds")

# Step 6: Calculate average delay per vehicle
# Average delay (d) is total delay divided by total vehicles per cycle.
d = D_total / N

print("\n--- Step 6: Final Average Delay ---")
print("The average deterministic delay per vehicle is calculated as:")
print(f"d = Total Delay / Total Vehicles = {D_total:.2f} veh-s / {N:.2f} veh")
print(f"Average delay per vehicle (d) = {d:.2f} seconds")

# Final Answer
final_answer = round(d, 2)
# print(f"\n<<< {final_answer} >>>")