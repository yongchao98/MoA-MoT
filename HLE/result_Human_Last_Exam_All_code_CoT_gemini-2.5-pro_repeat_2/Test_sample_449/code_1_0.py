import numpy as np

# Define the starting point
x0 = (3000, 4000)

# Calculate the Euclidean distance (norm) of the starting point from the origin
norm_x0 = np.sqrt(x0[0]**2 + x0[1]**2)

# Euler-Mascheroni constant
gamma = np.euler_gamma

# Calculate the potential kernel a(x0) using the asymptotic formula
# a(x) approx (2/pi)*log|x| + (2*gamma + 3*log(8))/pi
C = (2 * gamma + 3 * np.log(8)) / np.pi
a_x0 = (2 / np.pi) * np.log(norm_x0) + C

# The probability of hitting the set A is approximately 1/a(x0)
prob_hit = 1 / a_x0

# The probability of never hitting the set A is 1 - prob_hit
prob_escape = 1 - prob_hit

# Print the calculation steps
print(f"Starting point x0 = {x0}")
print(f"Distance from origin ||x0|| = {norm_x0}")
print(f"Value of the potential kernel a(x0) is approximately {a_x0:.4f}")
print(f"The probability of ever hitting the four neighbours is P(hit) ≈ 1/a(x0) = 1/{a_x0:.4f} ≈ {prob_hit:.4f}")
print("The probability of never hitting the four neighbours is P(escape) = 1 - P(hit)")
print(f"P(escape) ≈ 1 - {prob_hit:.4f} = {prob_escape:.4f}")

# Format the final answer to two significant digits
# A simple way to do this is to format the number in scientific notation and then reformat it.
# Or just use a format string that is reasonable.
final_answer = float(f"{prob_escape:.2g}")
print(f"\nThe approximate answer with two significant digits is: {final_answer}")
