import math

# Step 1: Determine u_1 from the properties of the welded sheet (B+C)
# Let 'a' be a symbolic variable for now.
# Sheet B is assumed to be a rectangle with vertices (0,0), (2a,0), (2a,a), (0,a) based on the problem statement.
# Area_B = 2a * a = 2a^2
# Mass_B = u_1 * Area_B = 2 * u_1 * a^2
# Centroid_B = (a, a/2)
#
# Sheet C is on top of B, same width 2a, height 4a.
# Area_C = 2a * 4a = 8a^2
# u_2 = 3
# Mass_C = u_2 * Area_C = 3 * 8a^2 = 24a^2
# Centroid_C = (a, a + 4a/2) = (a, 3a)
#
# The k-coordinate of the center of gravity, k_s, is given by:
# k_s = (k_B * Mass_B + k_C * Mass_C) / (Mass_B + Mass_C)
# We are given k_s = 2a.
# 2a = ((a/2)*(2*u_1*a^2) + (3a)*(24a^2)) / (2*u_1*a^2 + 24a^2)
# 2a = (u_1*a^3 + 72a^3) / (2*u_1*a^2 + 24a^2)
# 2a = a^3(u_1 + 72) / (2a^2(u_1 + 12))
# 2 = (u_1 + 72) / (2*(u_1 + 12))
# 4*(u_1 + 12) = u_1 + 72
# 4*u_1 + 48 = u_1 + 72
# 3*u_1 = 24
# u_1 = 8
u_1 = 8
print(f"Step 1: Calculated u_1 = {u_1}")

# Step 2: Determine 'a'
# The function f(x) and its derivatives are:
# f(x) = 0.5 * (ln(1 + x^4) + arctan(x^2))
# f'(x) = (2x^3 + x) / (1 + x^4)
# f''(x) = (-2x^6 - 3x^4 + 6x^2 + 1) / (1 + x^4)^2

x = 5

# Calculate f(5)
f5_val = 0.5 * (math.log(1 + x**4) + math.atan(x**2))
f5_rounded = round(f5_val, 1)

# Calculate f'(5)
f_prime5_val = (2*x**3 + x) / (1 + x**4)
f_prime5_rounded = round(f_prime5_val, 1)

# Calculate f''(5)
numerator = -2*x**6 - 3*x**4 + 6*x**2 + 1
denominator = (1 + x**4)**2
f_double_prime5_val = numerator / denominator
f_double_prime5_rounded = round(f_double_prime5_val, 1)

print(f"\nStep 2: Calculating f-terms and 'a'")
print(f"f(5) = {f5_val:.4f}, rounded to {f5_rounded}")
print(f"f'(5) = {f_prime5_val:.4f}, rounded to {f_prime5_rounded}")
print(f"f''(5) = {f_double_prime5_val:.4f}, rounded to {f_double_prime5_rounded}")

# Calculate the expression E = f(5) - 2f'(5) + 2f''(5)
E = f5_rounded - 2*f_prime5_rounded + 2*f_double_prime5_rounded
print(f"Expression E = {f5_rounded} - 2*({f_prime5_rounded}) + 2*({f_double_prime5_rounded}) = {E}")

# Calculate 'a'
a = (u_1 / 27) * (E**3)
print(f"a = ({u_1}/27) * ({E})^3 = {a}")

# Step 3: Determine relationship between 'l' and 'a'
# For sheet A, the y-coordinate of the center of gravity y_s must be 4a.
# The sheet is a trapezoid with vertices (0,0), (4a,0), (4a,4a), (0, 4a+l).
# The area is A = 16a^2 + 2al.
# The y_s is found by setting up the moment integral or by decomposition.
# From the derivation: 96a^2 + 12al = 48a^2 + 12al + l^2
# This simplifies to l^2 = 48a^2
print("\nStep 3: Deriving relationship between l and a")
print("For Sheet A, the condition y_s = 4a leads to the equation: l^2 = 48 * a^2")
print(f"So, l = sqrt(48) * a = 4 * sqrt(3) * a")

# Step 4: Calculate the final value of l
# Substitute the value of a into the expression for l.
l_val = 4 * math.sqrt(3) * a
print("\nStep 4: Final Calculation")
print(f"Substituting a = {a}:")
print(f"l = 4 * sqrt(3) * {a}")
print(f"l = {4*a} * sqrt(3)")
print(f"l = {l_val:.4f}")