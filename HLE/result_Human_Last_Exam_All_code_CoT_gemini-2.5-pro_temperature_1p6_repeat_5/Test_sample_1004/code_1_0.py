import sympy as sp

def solve_boat_problem():
    """
    Calculates the final downstream drift of the boat using symbolic mathematics.
    """
    # Define the symbolic variables
    L, v_0, v, y, k = sp.symbols('L v_0 v y k')

    # --- Step 1: Define the river flow velocity profile ---
    # The general form is v_flow = k * y * (L - y)
    v_flow_general = k * y * (L - y)

    # Use the condition v_flow(L/2) = v_0 to find k
    equation_for_k = sp.Eq(v_flow_general.subs(y, L/2), v_0)
    k_solution = sp.solve(equation_for_k, k)[0]

    # The specific flow velocity function
    v_flow = v_flow_general.subs(k, k_solution)

    # --- Step 2: Calculate the drift for the outbound trip (x_out) ---
    # The boat travels from y=0 to y=3*L/4.
    # The boat's velocity in the y-direction is v.
    # Therefore, dt = dy / v.
    # The drift in the x-direction is dx = v_flow * dt = (v_flow / v) * dy.
    # Integrate to find the total drift.
    integrand_out = v_flow / v
    x_out = sp.integrate(integrand_out, (y, 0, 3*L/4))

    # --- Step 3: Calculate the drift for the return trip (x_back) ---
    # The boat travels from y=3L/4 to y=0.
    # The drift calculation results in the same integral and value as the outbound trip.
    # dt = dy / v_y. Here v_y is -v. dy is negative. dt = dy / (-v).
    # dx = v_flow * dt = v_flow * (dy/-v)
    # Integrating from y=3L/4 to 0: ∫(v_flow * (dy/-v)) from 3L/4 to 0
    # = ∫(v_flow / v) from 0 to 3L/4. Same as x_out.
    x_back = x_out

    # --- Step 4: Calculate the total distance D ---
    total_distance = x_out + x_back

    # --- Step 5: Print the results clearly ---
    print("The velocity profile of the river is:")
    print(f"v_flow(y) = {v_flow}")
    print("-" * 40)
    
    print("The drift distance during the outbound trip (to y=3L/4) is:")
    print(f"x_out = {x_out}")
    print("-" * 40)

    print("The drift distance during the return trip (from y=3L/4) is the same:")
    print(f"x_back = {x_back}")
    print("-" * 40)

    print("The total distance between the start and return points is D = x_out + x_back.")
    
    # Get the coefficients and variables for the final expression
    num, den = total_distance.as_numer_denom()
    num_coeff = num.as_coeff_Mul()[0]
    num_vars = num.as_coeff_Mul()[1]
    den_coeff = den.as_coeff_Mul()[0]
    den_vars = den.as_coeff_Mul()[1]
    
    print("\nFinal Equation:")
    print(f"Distance = ({num_coeff} * {num_vars}) / ({den_coeff} * {den_vars})")

if __name__ == '__main__':
    solve_boat_problem()