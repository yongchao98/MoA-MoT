import numpy as np

# This script calculates the probability that a 2D simple random walk,
# conditioned to avoid the origin and starting from (3000, 4000),
# never enters the set of the origin's four neighbors.

# The formula for this probability is P = 1 - 2 / (ln|x0| + gamma + 1.5*ln(2)).

# Define parameters and constants
start_point = (3000, 4000)
gamma_const = np.euler_gamma

# Calculate the components of the formula
magnitude = np.sqrt(start_point[0]**2 + start_point[1]**2)
log_magnitude = np.log(magnitude)
ln2_term = 1.5 * np.log(2)

# Calculate the denominator for the hitting probability
denominator = log_magnitude + gamma_const + ln2_term

# Calculate the hitting and avoiding probabilities
p_hit = 2 / denominator
p_avoid = 1 - p_hit

# Output the equation and its step-by-step evaluation
print("The calculation for the probability P is as follows:")
print(f"P = 1 - 2 / (ln({magnitude}) + gamma + 1.5*ln(2))")
print(f"P \u2248 1 - 2 / ({log_magnitude:.3f} + {gamma_const:.3f} + {ln2_term:.3f})")
print(f"P \u2248 1 - 2 / ({denominator:.3f})")
print(f"P \u2248 1 - {p_hit:.3f}")
print(f"P \u2248 {p_avoid:.3f}")
print(f"\nThe final approximate answer with two significant digits is: {p_avoid:.2f}")