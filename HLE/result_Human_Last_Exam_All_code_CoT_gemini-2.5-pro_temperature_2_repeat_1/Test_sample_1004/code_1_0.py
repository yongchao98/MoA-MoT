import sympy

def calculate_boat_drift():
    """
    This function calculates the total downstream drift of a boat crossing a river,
    as described in the problem.
    It uses symbolic mathematics to derive the final expression.
    """
    # Step 1: Define symbolic variables for the problem.
    # L: width of the river
    # v_0: maximum flow velocity at the center
    # v: boat's speed relative to the water, perpendicular to the flow
    # y: variable for the distance from the starting bank
    L, v_0, v, y = sympy.symbols('L v_0 v y', real=True, positive=True)

    # Define the parabolic river flow velocity u(y).
    # The profile is u(y) = k*y*(L-y) where k is a constant.
    # We find k by using the condition u(L/2) = v_0.
    # v_0 = k * (L/2) * (L - L/2) = k * L**2 / 4  => k = 4*v_0/L**2
    river_velocity_u = (4 * v_0 / L**2) * y * (L - y)

    # The boat's cross-stream velocity component is v, so dy/dt = v.
    # This gives us a way to relate time and position: dt = dy / v.
    
    # Step 2: Calculate the drift for the first leg of the journey (outbound).
    # The boat travels from y=0 to y = L - L/4 = 3*L/4.
    # The drift dx is river_velocity_u * dt = river_velocity_u * (dy / v).
    # We integrate this expression from y=0 to y=3*L/4.
    drift_outbound = sympy.integrate(river_velocity_u / v, (y, 0, 3 * L / 4))

    # Step 3: Calculate the drift for the second leg of the journey (return).
    # The boat travels from y=3L/4 back to y=0.
    # The cross-stream velocity is now dy/dt = -v, so dt = -dy / v.
    # The drift dx is river_velocity_u * dt = river_velocity_u * (-dy / v).
    # We integrate from y=3*L/4 to y=0.
    drift_return = sympy.integrate(river_velocity_u * (-1 / v), (y, 3 * L / 4, 0))

    # Step 4: Calculate the total drift by summing the two parts.
    total_drift = drift_outbound + drift_return
    
    # Simplify the final expression.
    final_expression = sympy.simplify(total_drift)

    # Step 5: Output the result, highlighting the numbers in the final equation.
    print("The final expression for the total downstream distance is a fraction.")
    
    # We extract the numerator and denominator to find their numerical coefficients.
    numerator, denominator = final_expression.as_numer_denom()
    
    # as_coeff_mul() splits a term like 9*L*v_0 into its coefficient (9) and the rest (L*v_0).
    coeff_numerator = numerator.as_coeff_mul()[0]
    coeff_denominator = denominator.as_coeff_mul()[0]
    
    print(f"The number in the numerator is: {coeff_numerator}")
    print(f"The number in the denominator is: {coeff_denominator}")

    print("\nThe final equation for the total distance is:")
    sympy.pprint(final_expression)

# Execute the calculation and print the results.
calculate_boat_drift()