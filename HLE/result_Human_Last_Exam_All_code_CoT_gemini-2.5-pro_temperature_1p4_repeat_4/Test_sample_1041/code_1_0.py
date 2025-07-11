import math

# --- Parameters ---
# Starting position is (0, 300), so the distance from the origin is r.
r = 300.0
# The walk exits a disk of radius R.
R = 1000.0
# The target set is {(0,0), (2,0)}.
# The distance between the points is d.
d = 2.0
# The lattice spacing is a.
a = 1.0

# --- Calculation of Effective Radius ---
# The effective radius (r_eff) for two points is approximated by sqrt(a * d).
r_eff = math.sqrt(a * d)

# --- Calculation of the Probability ---
# The probability is given by the formula p = log(R/r) / log(R/r_eff).
# We use the natural logarithm (math.log).
numerator = math.log(R / r)
denominator = math.log(R / r_eff)
probability = numerator / denominator

# --- Output the result ---
print("This script calculates the probability that a 2D random walk starting at (0, 300)")
print("visits the set {(0,0), (2,0)} before leaving the disk of radius 1000.")
print("\nThe probability 'p' is approximated by the formula: p = ln(R/r) / ln(R/r_eff)")
print("\nWhere:")
print(f"  R = Radius of the exit disk = {R}")
print(f"  r = Starting distance from the origin = {r}")
print(f"  d = Distance between target points = {d}")
print(f"  a = Lattice spacing = {a}")
print(f"  r_eff = Effective radius of the target set = sqrt(a*d) = {r_eff:.3f}")

print("\nFinal equation:")
print(f"p = ln({R:.0f} / {r:.0f}) / ln({R:.0f} / {r_eff:.3f})")
print(f"p = ln({R/r:.3f}) / ln({R/r_eff:.3f})")
print(f"p = {numerator:.5f} / {denominator:.5f}")

print(f"\nThe calculated probability is: {probability:.3f}")

# Return the final answer in the required format
# print(f"\n<<<{probability:.3f}>>>")