import sympy as sp

# Define symbols
t, a = sp.symbols('t a', real=True, positive=True)

# Define the parametric equations
x = sp.cos(t)**3
y = sp.sin(t)**3

print("This script finds the value of 'a' for the arc defined by:")
print(f"x = {sp.pretty(x)}, y = {sp.pretty(y)}")
print(f"with the constraint 0 <= x <= a and an arc length of 3/2.\n")

# --- Step 1: Calculate derivatives for the arc length formula ---
dx_dt = sp.diff(x, t)
dy_dt = sp.diff(y, t)

print("Step 1: Calculate the derivatives dx/dt and dy/dt.")
print(f"dx/dt = {sp.pretty(dx_dt)}")
print(f"dy/dt = {sp.pretty(dy_dt)}\n")

# --- Step 2: Determine the integrand ---
integrand_sq = sp.simplify(dx_dt**2 + dy_dt**2)
integrand = sp.sqrt(integrand_sq)

print("Step 2: Simplify the expression for the integrand.")
print(f"(dx/dt)^2 + (dy/dt)^2 simplifies to: {sp.pretty(integrand_sq)}")
print(f"The integrand sqrt((dx/dt)^2 + (dy/dt)^2) is: {sp.pretty(integrand)}\n")

# --- Step 3: Determine the limits of integration ---
# The condition 0 <= x <= a for the arc in the first quadrant (where t is in [0, pi/2])
# corresponds to t being integrated from arccos(a^(1/3)) to pi/2.
# In this interval, sin(t) and cos(t) are non-negative, so the absolute value can be removed.
integrand_simplified = 3 * sp.sin(t) * sp.cos(t)
t_lower = sp.acos(a**(sp.S(1)/3))
t_upper = sp.pi/2

print("Step 3: Determine the integration limits and final integrand.")
print("For the arc in the first quadrant, the constraint 0 <= x <= a corresponds to")
print("t varying from t1 = arccos(a^(1/3)) to t2 = pi/2.")
print(f"In this region, the integrand simplifies to: {sp.pretty(integrand_simplified)}\n")

# --- Step 4: Evaluate the arc length integral ---
# The integral of 3*sin(t)*cos(t) is (3/2)*sin(t)^2.
# Evaluating from t_lower to t_upper:
# (3/2)*sin(pi/2)^2 - (3/2)*sin(acos(a^(1/3)))^2
# = (3/2) * 1 - (3/2) * (1 - a^(2/3))
# = (3/2) * a^(2/3)
arc_length_expr = sp.integrate(integrand_simplified, (t, t_lower, t_upper))

print("Step 4: Evaluate the arc length integral L(a).")
print(f"L(a) = Integral({sp.pretty(integrand_simplified)}, (t, {sp.pretty(t_lower)}, {sp.pretty(t_upper)}))")
print(f"The result of the integration is: L(a) = {sp.pretty(arc_length_expr)}\n")

# --- Step 5: Solve the equation L(a) = 3/2 ---
given_length = sp.S(3)/2
equation_to_solve = sp.Eq(arc_length_expr, given_length)

print("Step 5: Set the calculated arc length equal to the given length (3/2) and solve for 'a'.")
print("The final equation is:")

# We manually print the numbers in the equation as requested.
# The expression for arc length is (3/2) * a**(2/3). The given length is 3/2.
print("3/2 * a^(2/3) = 3/2")

solutions = sp.solve(equation_to_solve, a)

# Since 'a' is defined as a positive real number, we select the positive real solution.
final_a = solutions[0]

print(f"\nThe solution for 'a' is:")
print(f"a = {final_a}")
