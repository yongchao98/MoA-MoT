import sympy
from sympy import pi, Rational

# Announce the plan
print("The goal is to find the volume of the space enclosed by a cone S1 and an ellipsoid S2.")
print("The plan is as follows:")
print("1. Determine the plane where the cone is tangent to the ellipsoid.")
print("2. Calculate the volume of the cone (S1) from its vertex down to this plane of tangency.")
print("3. Calculate the volume of the cap of the ellipsoid (S2) above this plane of tangency.")
print("4. The enclosed volume is the difference between the volume of the cone part and the volume of the ellipsoid cap.")
print("-" * 20)

# Step 1: Find the plane of tangency
# The ellipsoid S2 is given by x^2/3 + y^2/4 + z^2/3 = 1.
# The cone S1 has its vertex at V = (0, 4, 0).
# The tangent plane to the ellipsoid at a point (x0, y0, z0) is x*x0/3 + y*y0/4 + z*z0/3 = 1.
# Since the vertex of the cone must lie on the tangent plane at the curve of tangency, we substitute its coordinates:
# 0*x0/3 + 4*y0/4 + 0*z0/3 = 1, which simplifies to y0 = 1.
y_tangency = 1
print(f"Step 1: The cone is tangent to the ellipsoid on the plane y = {y_tangency}.")
print("-" * 20)

# Step 2: Calculate the volume of the cone part (V_C)
# The base of the relevant cone part lies on the plane y=1.
# To find the radius of this base, we find the intersection of the ellipsoid and the plane y=1:
# x^2/3 + 1^2/4 + z^2/3 = 1  =>  x^2/3 + z^2/3 = 3/4  =>  x^2 + z^2 = 9/4.
# This is a circle with radius r = 3/2.
r_val = Rational(3, 2)
r_sq_val = r_val**2
# The height of this cone part is the distance from the vertex (y=4) to the base (y=1).
h_cone_val = 4 - 1
# The volume of a cone is V = (1/3)*pi*r^2*h.
V_C = (Rational(1, 3) * pi * r_sq_val * h_cone_val)
print(f"Step 2: The volume of the cone from its vertex (y=4) to its base at y=1 is calculated.")
print(f"   - Radius of the base r = {r_val}")
print(f"   - Height of the cone h = {h_cone_val}")
print(f"   - V_cone = (1/3) * \u03C0 * ({r_val})\u00B2 * {h_cone_val} = {V_C}")
print("-" * 20)


# Step 3: Calculate the volume of the ellipsoid cap (V_E)
# This is the volume of the part of the ellipsoid where y >= 1. The top of the ellipsoid is at y=2.
# We use the disk method. The area of a circular cross-section at a given y is A(y) = pi * R^2.
# From the ellipsoid equation, R^2 = x^2 + z^2 = 3 * (1 - y^2/4).
y = sympy.Symbol('y')
integrand = pi * 3 * (1 - y**2/4)
y_top_ellipsoid = 2
# We integrate this area from y=1 to y=2.
V_E = sympy.integrate(integrand, (y, y_tangency, y_top_ellipsoid))
print(f"Step 3: The volume of the ellipsoid cap for y >= 1 is calculated by integration.")
print(f"   - The integral is \u222B[{y_tangency},{y_top_ellipsoid}] (\u03C0 * 3 * (1 - y\u00B2/4)) dy")
print(f"   - V_ellipsoid_cap = {V_E}")
print("-" * 20)

# Step 4: Calculate the final volume
# The enclosed volume is V = V_C - V_E.
V_total = V_C - V_E
print("Step 4: The final enclosed volume is the difference between the two calculated volumes.")
print(f"Final Equation:")
print(f"Volume = V_cone - V_ellipsoid_cap")
print(f"Volume = ( (1/3) * \u03C0 * ({r_val})\u00B2 * {h_cone_val} ) - ( \u222B[{y_tangency},{y_top_ellipsoid}] (\u03C0 * 3 * (1 - y\u00B2/4)) dy )")
print(f"Volume = {V_C} - {V_E}")
print(f"Volume = {V_total}")