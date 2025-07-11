import math

# Problem parameters
R = 1000.0       # Radius of the stopping disk
S = (0.0, 300.0) # Starting point of the walk
T1 = (0.0, 0.0)  # First target point
T2 = (2.0, 0.0)  # Second target point
d_min = 1.0      # Lattice spacing, used as an effective radius for a point target

def distance(p1, p2):
    """Calculates Euclidean distance between two points."""
    return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)

# This problem is equivalent to finding the solution to a discrete Laplace equation
# with specific boundary conditions. We can approximate the solution using the
# corresponding continuous case (Brownian motion), where hitting probabilities
# are related to harmonic functions (potentials).

# The hitting probability u(z) for a target set {T1, T2} is approximated by:
# u(z) = A * (ln(R / |z - T1|)) + B * (ln(R / |z - T2|))
#
# To find the constants A and B, we enforce the boundary conditions u(T1)=1 and u(T2)=1.
# This leads to a system of two linear equations. Due to the symmetry of the problem, A = B.
# The resulting expression for A is:
# A = 1 / (ln(R / d_min) + ln(R / |T1 - T2|))

# Calculate required distances
dist_s_t1 = distance(S, T1)
dist_s_t2 = distance(S, T2)
dist_t1_t2 = distance(T1, T2)

# Step 1: Calculate the coefficient A
A_denominator = math.log(R / d_min) + math.log(R / dist_t1_t2)
A = 1.0 / A_denominator

# Step 2: Calculate the probability at the starting point S
term1 = math.log(R / dist_s_t1)
term2 = math.log(R / dist_s_t2)
probability = A * (term1 + term2)

# Step 3: Print the calculation step-by-step
print("The probability 'P' is calculated using an analytical approximation based on potential theory.")
print("The calculation proceeds as follows:")
print("-" * 60)
print("1. Define problem parameters:")
print(f"   Radius of disk, R = {R}")
print(f"   Starting point, S = {S}")
print(f"   Target set, T = {{{T1}, {T2}}}")
print("-" * 60)
print("2. The probability P at point S can be expressed as:")
print("   P = A * (ln(R / |S-T1|) + ln(R / |S-T2|))")
print("   where A is a coefficient determined by the boundary conditions.")
print("\n3. First, calculate the coefficient A:")
print("   A = 1 / (ln(R / 1) + ln(R / |T1-T2|))")
print(f"   A = 1 / (ln({R:.0f}) + ln({R:.0f} / {dist_t1_t2:.1f}))")
print(f"   A = 1 / ({math.log(R):.4f} + {math.log(R / dist_t1_t2):.4f})")
print(f"   A = 1 / ({A_denominator:.4f})")
print(f"   A = {A:.5f}")
print("-" * 60)
print("4. Now, substitute all values into the probability equation:")
print(f"   P = {A:.5f} * (ln({R:.0f} / {dist_s_t1:.1f}) + ln({R:.0f} / {dist_s_t2:.3f}))")
print(f"   P = {A:.5f} * ({term1:.4f} + {term2:.4f})")
print(f"   P = {A:.5f} * ({(term1 + term2):.4f})")
print(f"   P = {probability:.5f}")
print("-" * 60)
print(f"The final probability, rounded to three significant digits, is {probability:.3g}.")
