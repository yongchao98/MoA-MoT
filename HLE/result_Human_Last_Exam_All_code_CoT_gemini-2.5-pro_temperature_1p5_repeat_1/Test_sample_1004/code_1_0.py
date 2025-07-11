import sympy

# 1. Define the symbolic variables
# L: width of the river
# v_0: maximum flow velocity at the center
# v: boat's speed relative to the water
# y: distance from the starting bank
L, v_0, v, y = sympy.symbols('L v_0 v y', positive=True, real=True)

# 2. Define the river's flow velocity as a piecewise function
# The velocity profile is triangular: 0 at the banks (y=0, y=L), and v_0 at the center (y=L/2).
v_flow_part1 = (2 * v_0 / L) * y
v_flow_part2 = (2 * v_0 / L) * (L - y)
v_flow = sympy.Piecewise(
    (v_flow_part1, y <= L/2),
    (v_flow_part2, y > L/2)
)

# 3. Calculate the downstream drift
# The integrand for drift is derived from dx = v_flow(y) * dt, where dt = dy / v.
# So, we integrate (v_flow(y) / v) with respect to y.

# The boat travels from y=0 to y = L - L/4 = 3L/4 on the outbound trip.
drift_outbound = sympy.integrate(v_flow / v, (y, 0, 3*L/4))

# The drift on the return trip is identical because the boat traverses the same y-path
# and the flow velocity only depends on y.
drift_inbound = drift_outbound

# The total distance is the sum of the drift from both trips.
total_distance = drift_outbound + drift_inbound

# 4. Simplify the final expression to get the formula
total_distance_simplified = sympy.simplify(total_distance)

# 5. Extract the numbers from the simplified symbolic expression for clear output.
# The result is of the form (numerator/denominator) * (L*v_0/v).
if total_distance_simplified.is_Mul:
    coefficient, _ = total_distance_simplified.as_coeff_mul()
    numerator = sympy.numer(coefficient)
    denominator = sympy.denom(coefficient)
else:
    # Handle the case where the expression is simpler than expected
    numerator = total_distance_simplified
    denominator = 1

print("The problem asks for the total downstream distance the boat travels.")
print("This is the sum of the downstream drift on the trip out and the trip back.")
print(f"\nDrift on the outbound trip (from y=0 to y=3L/4): {drift_outbound}")
print(f"Drift on the return trip (from y=3L/4 to y=0): {drift_inbound}")
print(f"\nTotal drift = (Outbound Drift) + (Inbound Drift)")
print(f"Total drift = {total_distance_simplified}\n")

print("Final Answer:")
print("The distance between the boat's returning position on the home bank and its original starting point is:")
# Output the final equation with the calculated numbers
print(f"Distance = ({numerator} * L * v_0) / ({denominator} * v)")