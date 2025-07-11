import math

# Define the constants
# Euler-Mascheroni constant
gamma_E = 0.5772156649
# Natural log of sqrt(8)
log_sqrt8 = math.log(8) / 2
# Combined constant
C = gamma_E + log_sqrt8

# Define the start and target points
x_start = (3000, 4000)
# A is the set of neighbours of the origin, which are at distance 1
r_A = 1
# Calculate the distance of the starting point from the origin
r_start = math.sqrt(x_start[0]**2 + x_start[1]**2)

# Define the potential kernel function h(r) for large r,
# normalized such that h(1) = -1.
# h(r) is approximately -2 * (ln(r) + C)
# The value at r=1 for the normalized kernel is exactly -1 by definition.
h_A = -1

# Calculate the value of h at the starting point using the asymptotic formula.
h_start = -2 * (math.log(r_start) + C)

# The probability of ever hitting the set A is h_A / h_start.
prob_hit_A = h_A / h_start

# The probability of never hitting A is 1 minus the probability of hitting A.
prob_never_hit_A = 1 - prob_hit_A

# Print the calculation steps
print(f"The starting position is x_0 = {x_start}, with distance |x_0| = {r_start:.0f} from the origin.")
print(f"The set to avoid is A, the neighbors of the origin, with distance |y| = {r_A} for y in A.")
print("The probability of never hitting A is given by the formula P = 1 - h(A)/h(x_0).")
print(f"With our chosen normalization for the potential kernel h(x), we have h(A) = {h_A}.")
print(f"The value of h(x_0) is calculated using its asymptotic form: h({r_start:.0f}) ≈ -2 * (ln({r_start:.0f}) + γ_E + ln(√8))")
print(f"h({r_start:.0f}) ≈ -2 * ({math.log(r_start):.4f} + {gamma_E:.4f} + {log_sqrt8:.4f}) ≈ {h_start:.4f}")
print(f"P = 1 - ({h_A} / {h_start:.4f})")
print(f"P ≈ {prob_never_hit_A:.4f}")

# Final answer with two significant digits
print(f"The approximate probability is {prob_never_hit_A:.2f}.")
<<<0.95>>>