import math

# Step 1: Define Given Parameters
s_hr = 2160  # Saturation flow rate (veh/hour)
v_hr = 600   # Approach average flow rate (veh/hour)
R_disp = 56  # Displayed red time (seconds)
Y = 3        # Displayed yellow time (seconds)
AR = 2       # All-red time (seconds)
g = 30       # Effective green time (seconds)
t_L = 4      # Total lost time (seconds)
percent_green_arrival = 0.40 # 40% of traffic arrives during green
percent_red_arrival = 0.60   # 60% of traffic arrives during red

# Step 2: Convert Units
# Convert flow rates from veh/hour to veh/sec
s = s_hr / 3600
v = v_hr / 3600

# Step 3: Calculate Cycle Components
# Effective red time (r) is the sum of displayed red and all-red time
r = R_disp + AR
# Cycle length (C) is the sum of effective green and effective red time
C = g + r

# Step 4: Calculate Arrival Characteristics
# Total vehicles arriving per cycle (N)
N = v * C
# Number of vehicles arriving during the red interval (N_r)
N_r = N * percent_red_arrival
# The queue at the end of red (Q_r) is equal to the number of arrivals during red
Q_r = N_r
# Number of vehicles arriving during the green interval (N_g)
N_g = N * percent_green_arrival
# Arrival rate during the green interval (lambda1)
lambda1 = N_g / g

# Step 5: Analyze the Queue
# Time required for the queue to clear during green (t_c)
# This assumes the intersection is undersaturated (s > lambda1)
if s <= lambda1:
    print("Error: Saturation flow rate is not greater than arrival rate during green. Intersection is oversaturated.")
else:
    t_c = Q_r / (s - lambda1)

    # Step 6: Calculate Delay
    # Total delay per cycle (D_total) using the area of the delay polygon
    total_delay = 0.5 * Q_r * (r + t_c)

    # Average delay per vehicle (d)
    avg_delay = total_delay / N

    # Step 7: Format and Print the Output
    print("The average deterministic delay per vehicle is calculated as follows:")
    print(f"Average Delay = (0.5 * Queue_at_end_of_red * (Effective_red_time + Queue_clearance_time)) / Total_vehicles_per_cycle")
    print(f"Average Delay = (0.5 * {Q_r:.2f} veh * ({r:.2f} s + {t_c:.2f} s)) / {N:.2f} veh")
    print(f"Average Delay = {avg_delay:.2f} seconds")

    # Final Answer
    print(f"\n<<< {avg_delay:.2f} >>>")
