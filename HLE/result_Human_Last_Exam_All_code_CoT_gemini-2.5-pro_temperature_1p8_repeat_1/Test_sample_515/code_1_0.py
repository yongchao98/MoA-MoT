import math

# The problem asks for the probability that a specific random walk starting
# from (0,1) hits the set of neighbors of (3600,0).
# The walk is a 2D simple random walk conditioned to never visit the origin.
# This conditioned walk is a Doob's h-transform process with respect to
# the potential kernel h(x) = a(x) of the 2D SRW.

# The probability of this process starting from x_0 hitting a distant set A
# (approximated by a point z) is given by a(x_0) / a(z).

# 1. The starting point is x_0 = (0,1).
# The potential kernel value a(x_0) = a(0,1) is exactly 1 by standard results.
a_x0 = 1.0

# 2. The target set is the neighbors of z = (3600,0).
# We approximate the potential kernel value for this set by a(z).
# For large distances ||x||, a(x) has a known asymptotic formula:
# a(x) ~= (2/pi) * (ln(||x||) + gamma + log(8))
# where gamma is the Euler-Mascheroni constant.

# Define constants
gamma = math.euler
log_8 = math.log(8)
z_norm = 3600.0

# Calculate the term in the parenthesis
log_z = math.log(z_norm)
term_in_parenthesis = log_z + gamma + log_8

# Calculate a(z)
a_z = (2 / math.pi) * term_in_parenthesis

# The probability is the ratio a(x_0) / a(z)
probability = a_x0 / a_z

# We need to output the approximate answer with two significant digits.
# The calculation shows:
# log(3600) ~= 8.188786
# gamma ~= 0.577216
# log(8) ~= 2.079441
# sum = 10.845443
# a(z) = (2/pi) * sum ~= 6.9052
# probability = 1 / 6.9052 ~= 0.1448
# Rounding to two significant digits gives 0.14.

print("This plan calculates the probability that a 2D random walk conditioned to avoid the origin, starting at (0,1), reaches the neighborhood of (3600,0).")
print("The probability 'P' can be approximated by the ratio of the potential kernel 'a(x)' at the start and target points.")
print("P ≈ a(x_0) / a(z)")
print("where x_0 = (0,1) and z = (3600,0).")
print("\nStep 1: Determine a(x_0).")
print("The starting point is x_0 = (0,1). From standard results for 2D SRW, a(0,1) = 1.")
print(f"a(x_0) = {a_x0}")

print("\nStep 2: Calculate a(z) for z = (3600,0).")
print("For a point x far from the origin, a(x) is approximated by: a(x) ≈ (2/π) * (ln(||x||) + γ + ln(8))")
print(f"||z|| = {z_norm}")
print(f"ln(||z||) = ln({z_norm}) ≈ {log_z}")
print(f"γ (Euler's constant) ≈ {gamma}")
print(f"ln(8) ≈ {log_8}")
print(f"The sum in the parenthesis is ≈ {term_in_parenthesis}")
print(f"a(z) ≈ (2/π) * {term_in_parenthesis} ≈ {a_z}")

print("\nStep 3: Calculate the final probability.")
print(f"P ≈ a(x_0) / a(z) = {a_x0} / {a_z} ≈ {probability}")

# The question asks for the answer with two significant digits.
# The hundredths place is the second significant digit.
final_answer = round(probability, 2)
print(f"\nThe approximate probability is {probability:.4f}, which is {final_answer} to two significant digits.")

# Final explicit statement for the answer format.
print("\nFinal calculation:")
print(f"Probability = {a_x0} / {a_z}")