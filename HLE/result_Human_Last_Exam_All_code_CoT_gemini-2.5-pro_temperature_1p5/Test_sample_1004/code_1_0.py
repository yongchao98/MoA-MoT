import sympy

def solve_river_crossing_problem():
    """
    This function symbolically solves the river boat problem to find the total
    downstream displacement. It uses the sympy library to perform symbolic integration.
    """
    # Define the symbolic variables for the problem.
    # L: width of the river
    # v0_sym: maximum flow velocity at the center. We use a symbol with a subscript for display.
    # v: boat's speed relative to the water (perpendicular to the flow)
    # y: distance from the starting bank
    L, v = sympy.symbols('L v', positive=True, real=True)
    v0_sym = sympy.Symbol('v_0', positive=True, real=True)
    y = sympy.symbols('y', real=True)

    # --- Step 1: Define the river's flow velocity profile u(y) ---
    # The profile is parabolic, zero at y=0 and y=L, and v0 at y=L/2.
    # We can model this as u(y) = k * y * (L - y).
    # To find the constant k, we use the condition that u(L/2) = v_0:
    # v_0 = k * (L/2) * (L - L/2) = k * L**2 / 4
    # Solving for k gives: k = 4 * v_0 / L**2
    k = 4 * v0_sym / L**2
    u_y = k * y * (L - y)
    print("Step 1: The river flow velocity profile is derived as:")
    print(f"u(y) = {u_y}\n")

    # --- Step 2: Calculate the downstream drift during the outward journey ---
    # The boat travels from y=0 to the turning point y = L - L/4 = 3*L/4.
    # The boat's velocity across the river is dy/dt = v, which means dt = dy/v.
    # The downstream drift is the integral of the flow velocity over time.
    # dx_outward = integral(u(y) * dt) = integral(u(y)/v * dy) from y=0 to y=3L/4.
    y_turn = sympy.Rational(3, 4) * L
    integrand = u_y / v
    drift_outward = sympy.integrate(integrand, (y, 0, y_turn))
    print("Step 2: The downstream drift during the outward journey (from y=0 to y=3L/4) is:")
    print(f"Drift_outward = {drift_outward}\n")

    # --- Step 3: Calculate the downstream drift during the return journey ---
    # The boat travels from y=3L/4 back to y=0.
    # The boat's velocity across the river is dy/dt = -v, so dt = dy/(-v).
    # The downstream drift dx_return = integral(u(y) * dt) from y=3L/4 to y=0.
    # This integral is the same value as the outward drift.
    drift_return = sympy.integrate(u_y / (-v), (y, y_turn, 0))
    print("Step 3: The downstream drift during the return journey (from y=3L/4 to y=0) is:")
    print(f"Drift_return = {drift_return}\n")
    
    # --- Step 4: Calculate the total displacement ---
    # The total displacement is the sum of the outward and return drifts.
    total_displacement = drift_outward + drift_return
    print("Step 4: The total displacement is the sum of the outward and return drifts.\n")
    
    # --- Step 5: Format and print the final result ---
    # The symbolic result from sympy is of the form (numerator/denominator) * (L*v_0/v).
    # We parse the expression to extract the numerical coefficients.
    coeff = total_displacement.as_coeff_Mul()[0]
    
    # The coefficient will be a sympy Fraction object, e.g., sympy.Rational(9, 8)
    num = coeff.p
    den = coeff.q
    
    print("--- FINAL RESULT ---")
    print("The total distance between the boat's returning position and its original starting point is given by the following equation:")
    # We construct the final equation string as requested, showing each number.
    print(f"Total Displacement = ({num} * L * v_0) / ({den} * v)")

# Execute the function to solve the problem and print the results.
solve_river_crossing_problem()