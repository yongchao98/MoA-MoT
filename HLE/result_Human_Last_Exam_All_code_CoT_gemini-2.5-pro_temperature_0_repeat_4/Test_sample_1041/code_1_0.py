import math

# Define the parameters of the problem
R = 1000.0  # Radius of the disk
start_point = (0, 300)
target_points = [(0, 0), (2, 0)]

# 1. Determine the center of the target set
center_x = (target_points[0][0] + target_points[1][0]) / 2
center_y = (target_points[0][1] + target_points[1][1]) / 2
target_center = (center_x, center_y)

# 2. Calculate the distance from the start point to the center of the target set
r_s_sq = (start_point[0] - target_center[0])**2 + (start_point[1] - target_center[1])**2
r_s = math.sqrt(r_s_sq)

# 3. Calculate the effective radius of the target set
dist_sq = (target_points[0][0] - target_points[1][0])**2 + (target_points[0][1] - target_points[1][1])**2
dist = math.sqrt(dist_sq)
epsilon_eff = math.sqrt(dist)

# 4. Apply the harmonic measure formula
# P = log(R / r_s) / log(R / epsilon_eff)
# This can be rewritten using log properties as:
# P = (log(R) - log(r_s)) / (log(R) - log(epsilon_eff))
# P = (log(R) - 0.5*log(r_s^2)) / (log(R) - 0.5*log(dist))

log_R = math.log(R)
log_r_s = math.log(r_s)
log_epsilon_eff = math.log(epsilon_eff)

probability = (log_R - log_r_s) / (log_R - log_epsilon_eff)

# Print the components of the calculation
print("Problem Parameters:")
print(f"Radius of the disk, R = {R}")
print(f"Starting point, S = {start_point}")
print(f"Target set, A = {target_points}")
print("-" * 30)
print("Derived Values for the Formula:")
print(f"Center of target set, a_c = {target_center}")
print(f"Distance from start to target center, r_s = sqrt({r_s_sq}) = {r_s:.4f}")
print(f"Distance between target points, d = {dist}")
print(f"Effective radius of target set, epsilon_eff = sqrt(d) = {epsilon_eff:.4f}")
print("-" * 30)
print("Final Equation:")
print(f"P = log(R / r_s) / log(R / epsilon_eff)")
print(f"P = (log({R}) - log({r_s:.4f})) / (log({R}) - log({epsilon_eff:.4f}))")
print(f"P = ({log_R:.4f} - {log_r_s:.4f}) / ({log_R:.4f} - {log_epsilon_eff:.4f})")
numerator = log_R - log_r_s
denominator = log_R - log_epsilon_eff
print(f"P = {numerator:.4f} / {denominator:.4f}")
print("-" * 30)
print(f"The probability is: {probability:.3f}")
