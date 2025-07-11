import sympy as sp

# Define symbols for symbolic mathematics
y = sp.Symbol('y')
pi = sp.pi

# Step 1 & 2: From the geometric analysis, the plane of tangency is y=1.
# The common base is the circle x^2 + z^2 = 9/4 at y=1.
# The radius squared of this base is r_sq.
r_sq = sp.S(9)/4

# Step 3: Calculate the volume of the cone segment (V_cone).
# The cone segment has its vertex at y=4 and base at y=1.
# Height of the cone, h_cone = 4 - 1 = 3.
h_cone = 3
# Volume of a cone is (1/3) * pi * r^2 * h.
V_cone = (sp.S(1)/3) * pi * r_sq * h_cone

# Step 4: Calculate the volume of the ellipsoid cap (V_cap).
# This is found by integrating the cross-sectional area A(y) from y=1 to y=2.
# The cross-sectional area of the ellipsoid at height y is A(y) = pi * 3 * (1 - y^2/4).
Area_y = pi * (sp.S(3)/4) * (4 - y**2)
# Integrate A(y) from y=1 to y=2 (the top of the ellipsoid).
V_cap = sp.integrate(Area_y, (y, 1, 2))

# Step 5: The final enclosed volume is the difference between the two volumes.
V_final = V_cone - V_cap

# --- Output the Final Equation ---
# The problem asks to output each number in the final equation.
# The equation is V_final = V_cone - V_cap.
print("The final volume is calculated by the equation: V = V_cone - V_ellipsoid_cap")
print("\nSubstituting the calculated values:")
print(f"V = {V_cone} - {V_cap} = {V_final}")

# For a final numerical answer
# print(f"\nThe numerical value of the enclosed volume is: {V_final.evalf()}")
