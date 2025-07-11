import sympy

def solve_river_crossing():
    """
    Calculates the total downstream drift of a boat crossing a river,
    traveling partway across, and then returning.
    """
    # Step 1: Define the symbolic variables
    L, v0, v, y = sympy.symbols('L v_0 v y', real=True, positive=True)

    # Step 2: Define the river flow velocity as a function of distance y from the shore.
    # The flow is v_f(y) = k*y*(L-y).
    # At the center, y = L/2, the velocity is v0.
    # v0 = k * (L/2) * (L - L/2) = k * L**2 / 4  =>  k = 4*v0 / L**2
    v_f = (4 * v0 / L**2) * y * (L - y)

    # Step 3: Calculate the drift on the way out.
    # The boat travels from y=0 to y=3L/4.
    # The boat's speed perpendicular to the current is v.
    # dt = dy/v, so dx = v_f * dt = (v_f/v) * dy.
    integrand_out = v_f / v
    drift_out = sympy.integrate(integrand_out, (y, 0, 3 * L / 4))

    # Step 4: Calculate the drift on the way back.
    # The boat travels from y=3L/4 back to y=0.
    # The boat's speed perpendicular to the current is -v.
    # dt = dy/(-v), so dx = v_f * dt = (v_f/(-v)) * dy.
    integrand_back = v_f / (-v)
    drift_back = sympy.integrate(integrand_back, (y, 3 * L / 4, 0))

    # Step 5: Calculate the total drift by summing the two parts.
    total_drift = drift_out + drift_back

    # Step 6: Extract the numerical coefficient from the final simplified expression.
    # The expression will be of the form: C * L * v0 / v
    # We can isolate the numerical coefficient C.
    numerical_coeff = total_drift / (L * v0 / v)
    
    # Get the numerator and denominator of the coefficient
    numerator = sympy.numer(numerical_coeff)
    denominator = sympy.denom(numerical_coeff)

    # Step 7: Print the final results.
    print("The total distance between the returning position and the original starting point is given by the equation:")
    # Use sympy.pretty_print for a nicer output format if available, otherwise just print
    try:
        sympy.pretty_print(sympy.Eq(sympy.Symbol('Total_Drift'), total_drift))
    except ImportError:
        print(f"Total_Drift = {total_drift}")

    print("\nThe final simplified equation for the total drift is:")
    print(f"Total_Drift = ({numerator} * L * v_0) / ({denominator} * v)")
    
    print("\nThe numbers in the final equation are:")
    print(f"Numerator coefficient: {numerator}")
    print(f"Denominator coefficient: {denominator}")

solve_river_crossing()
<<<9*L*v_0/(8*v)>>>