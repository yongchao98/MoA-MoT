import sympy

# Define the symbolic variables for river width (L), max flow velocity (v_0),
# and the boat's relative speed (v). y is the integration variable.
L, v0, v, y = sympy.symbols('L v0 v y', positive=True, real=True)

# 1. Model the river's flow velocity v_x(y)
# The flow is zero at the banks (y=0, y=L) and v0 at the center (y=L/2).
# It is proportional to the distance from the nearest shore.
# This describes a triangular (piecewise linear) velocity profile.
v_x = sympy.Piecewise(
    ((2 * v0 / L) * y, y <= L / 2),
    ((2 * v0 / L) * (L - y), y > L / 2)
)

# 2. Set up the integral for the downstream displacement (drift)
# The boat's speed across the river is a constant v.
# The time to travel a vertical distance dy is dt = dy / v.
# The drift dx during that time is v_x(y) * dt.
# The integrand for the total drift is therefore v_x(y) / v.
integrand = v_x / v

# 3. Calculate the drift for the outward and return trips.
# Outward trip: The boat travels from the home bank (y=0) to a point
# L/4 from the opposite bank, which is y = L - L/4 = 3*L/4.
dx_outward = sympy.integrate(integrand, (y, 0, 3 * L / 4))

# Return trip: The boat travels from y = 3*L/4 back to y = 0.
# The magnitude of the drift over the same y-interval is identical.
dx_return = dx_outward

# 4. Calculate the total displacement.
# It's the sum of the drift from the outward and return trips.
total_displacement = dx_outward + dx_return

# Simplify the final expression
simplified_displacement = sympy.simplify(total_displacement)

# Extract the numerical coefficient from the final formula for clarity.
# The formula will be in the form: Coefficient * (L * v0 / v)
coefficient = simplified_displacement / (L * v0 / v)
numerator, denominator = sympy.fraction(coefficient)

# Print the final results
print("The final formula for the total downstream distance is derived from integrating the boat's drift.")
print("\nCalculating...\n")
print("The total distance between the start and end points is:")
# The unicode character for v_0 is u"\u_0".
print(f"Total Distance = ({numerator}/{denominator}) * (L * v\u2080 / v)")

print("\nThe numbers in the final equation are:")
print(f"Numerator: {numerator}")
print(f"Denominator: {denominator}")