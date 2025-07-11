import math

# Parameters of the problem
R = 1000  # Radius of the disk
x_start, y_start = 0, 300  # Starting point S
target_set = [(0, 0), (2, 0)]  # Target set A

# The center of the target set A is (1, 0)
x_center, y_center = 1, 0

# The distance 'r' from the starting point to the center of the target set
r = math.sqrt((x_start - x_center)**2 + (y_start - y_center)**2)

# The distance 'd' between the two points in the target set
d = math.sqrt((target_set[0][0] - target_set[1][0])**2 + (target_set[0][1] - target_set[1][1])**2)

# Effective radius 'epsilon' of a single lattice point for SRW on Z^2
epsilon = 0.25

# The effective radius 'rho' of the two-point target set
rho = math.sqrt(d * epsilon)

# The probability is given by the formula P = ln(R/r) / ln(R/rho)
# which is equivalent to (ln(R) - ln(r)) / (ln(R) - ln(rho))
log_R = math.log(R)
log_r = math.log(r)
log_rho = math.log(rho)

probability = (log_R - log_r) / (log_R - log_rho)

# Output the steps of the calculation
print("Calculation Steps:")
print(f"1. Outer radius R = {R}")
print(f"2. Starting distance r = |S - center(A)| = |({x_start},{y_start}) - ({x_center},{y_center})| = {r:.4f}")
print(f"3. Distance between target points d = |(0,0) - (2,0)| = {d}")
print(f"4. Effective radius of a single point epsilon = {epsilon}")
print(f"5. Effective radius of the target set rho = sqrt(d * epsilon) = sqrt({d} * {epsilon}) = {rho:.4f}")
print("\nThe probability P is calculated using the formula: (ln(R) - ln(r)) / (ln(R) - ln(rho))")
print(f"P = (ln({R}) - ln({r:.4f})) / (ln({R}) - ln({rho:.4f}))")
print(f"P = ({log_R:.4f} - {log_r:.4f}) / ({log_R:.4f} - {log_rho:.4f})")
print(f"P = {log_R - log_r:.4f} / {log_R - log_rho:.4f}")
print(f"P = {probability:.4f}")

# Final answer with three significant digits
print(f"\nThe probability, rounded to three significant digits, is: {probability:.3g}")
