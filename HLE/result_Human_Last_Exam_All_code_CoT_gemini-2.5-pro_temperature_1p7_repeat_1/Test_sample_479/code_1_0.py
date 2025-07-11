#
# Step-by-step derivation for the maximal genus of the surface.
#

# Let S be the surface, S = ∂Ω. Let g be the genus of S.
# Let H be the mean curvature of S with respect to the outward normal.

# Condition 1: S is the boundary of a compact region in R^3.
# This implies S is a compact, connected, orientable surface.

# Condition 2: The mean curvature H is never zero on S.
# Since S is connected and H is continuous, H must be either always positive or always negative.

# Theorem 1: Any compact surface in R^3 (like S) must contain at least one point p
# where H(p) > 0. This is shown by considering the point on the surface furthest
# from the origin and comparing it with the enclosing sphere tangent at that point.
# The mean curvature of the surface must be at least that of the sphere (which is positive).
# So, there is a point p with H(p) > 0.

# Conclusion 1: Combining Condition 2 and Theorem 1, we must have H > 0 everywhere on S.
positive_mean_curvature_everywhere = True

# Theorem 2 (Meeks-Simon-Yau, Schoen-Yau): A compact, embedded surface in R^3
# with H > 0 everywhere must be topologically a sphere.
must_be_a_sphere = positive_mean_curvature_everywhere

# Property: The genus of a sphere is 0.
if must_be_a_sphere:
    genus = 0
else:
    # This case is not possible under the problem's conditions.
    genus = "undefined"

# The only possible genus is 0. Therefore, the maximal genus is 0.
maximal_genus = 0

# Final Answer Equation
# The variables are for explanatory purposes, the core logic is in the theorems.
g_max = maximal_genus
print(f"Let g be the genus of the surface.")
print(f"From the given conditions, the mean curvature H must be positive everywhere.")
print(f"A surface with H > 0 everywhere must be a sphere, which has genus 0.")
print(f"Therefore, the only possible genus is 0.")
print(f"The maximal genus is g_max = {g_max}")
