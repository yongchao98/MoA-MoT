import sympy

# Define the symbolic variables for the problem.
# L: width of the river
# v_0: maximum flow velocity at the center
# v: boat's speed relative to the water
# y: distance from the starting bank
L, v_0, v, y = sympy.symbols('L v_0 v y', positive=True, real=True)

# 1. Define the river's flow velocity profile as a function of y.
# The profile is parabolic: v_flow = k*y*(L-y).
# We find k from the condition v_flow(L/2) = v_0, which gives k = 4*v_0/L**2.
v_flow = (4 * v_0 / L**2) * y * (L - y)

# 2. Calculate the drift for the first leg of the journey (outward).
# The boat travels from y=0 to y=3L/4 with a perpendicular velocity v.
# The time to cross dy is dt = dy/v.
# The drift dx during this time is v_flow * dt.
# So, the integrand for the drift is v_flow / v.
drift_integrand_leg1 = v_flow / v
drift_leg1 = sympy.integrate(drift_integrand_leg1, (y, 0, 3 * L / 4))

# 3. Calculate the drift for the second leg of the journey (return).
# The boat travels from y=3L/4 back to y=0 with a perpendicular velocity -v.
# The time to cross dy is dt = dy/(-v).
# The integrand is v_flow / (-v). We integrate from the starting y to the final y.
drift_integrand_leg2 = v_flow / (-v)
drift_leg2 = sympy.integrate(drift_integrand_leg2, (y, 3 * L / 4, 0))

# 4. The total distance is the sum of the drifts from both legs.
total_distance = drift_leg1 + drift_leg2

# 5. The problem asks to output the final equation with its numbers.
# We will extract the numerical coefficients from the symbolic result.
# The result is of the form (numerator_coefficient * variables) / (denominator_coefficient * variables)
num, den = total_distance.as_numer_denom()

# Extract the numerical coefficient from the numerator part.
num_coeff = num.as_coeff_mul(L, v_0)[0]
# Extract the numerical coefficient from the denominator part.
den_coeff = den.as_coeff_mul(v)[0]

# 6. Print the final result in the requested format.
print("The total distance between the returning position and the starting point is given by the equation:")
print(f"Total Distance = ({num_coeff} * L * v_0) / ({den_coeff} * v)")