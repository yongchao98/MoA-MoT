import sympy as sp

def solve_river_crossing():
    """
    Calculates the total downstream drift of a boat crossing a river with a parabolic flow profile.
    """
    # Define the symbolic variables for the problem
    L, v0, v, y, k = sp.symbols('L v_0 v y k')

    # Step 1 & 2: Determine the river's velocity profile u(y)
    # The general form is u(y) = k * y * (L - y)
    # We find k using the condition that velocity is v0 at the center y = L/2
    equation_for_k = sp.Eq(k * (L/2) * (L - L/2), v0)
    k_solution = sp.solve(equation_for_k, k)[0]
    
    # The specific velocity profile of the river
    u_y = k_solution * y * (L - y)

    # Step 3: Calculate the downstream drift for the outbound trip (x1)
    # The boat travels from y=0 to y=3L/4.
    # The time element is dt = dy / v. The drift is dx = u(y) * dt.
    # We integrate dx to find the total drift x1.
    integrand_outbound = u_y / v
    x1 = sp.integrate(integrand_outbound, (y, 0, 3*L/4))

    # Step 4: Calculate the downstream drift for the inbound trip (x2)
    # The boat travels from y=3L/4 back to y=0.
    # The time element is dt = dy / (-v). The drift is dx = u(y) * dt.
    # The integral is from y=3L/4 to y=0.
    # Note: The integral from a to b of -f(x) is the same as the integral from b to a of -f(x),
    # which is the same as the integral from a to b of f(x). So x2 = x1.
    x2 = x1

    # Step 5: Calculate the total downstream distance
    total_drift = x1 + x2

    # Step 6: Print the final result clearly
    print("The total distance between the boat's returning position and its original starting point is given by the formula:")
    
    # Extract the numerical coefficients from the final expression
    num, den = total_drift.as_numer_denom()
    num_coeff = num.as_coeff_mul()[0]
    den_coeff = den.as_coeff_mul()[0]
    
    # Print the final equation showing each number
    print(f"\nDistance = ({num_coeff} * v_0 * L) / ({den_coeff} * v)")

solve_river_crossing()