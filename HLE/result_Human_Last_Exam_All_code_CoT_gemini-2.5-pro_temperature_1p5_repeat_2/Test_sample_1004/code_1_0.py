import sympy

def solve_boat_problem():
    """
    Calculates the total downstream displacement of a boat crossing a river.
    """
    # Step 1: Define the symbolic variables and the river velocity profile.
    # L: width of the river
    # v_0: maximum velocity of the river at the center
    # v: speed of the boat relative to the water, perpendicular to the flow
    # y: distance from the starting bank
    L, v_0, v, y, k = sympy.symbols('L v_0 v y k')

    print("Step 1: Determine the river's velocity profile v_r(y).")
    # The velocity is zero at y=0 and y=L, and proportional to distance from shores.
    # A parabolic profile fits this: v_r(y) = k * y * (L - y)
    v_r_general = k * y * (L - y)
    print(f"Assuming a parabolic velocity profile: v_r(y) = k * y * (L - y)")

    # We are given that the velocity is v_0 at the center (y = L/2).
    # We use this condition to solve for the constant k.
    equation_for_k = sympy.Eq(v_r_general.subs(y, L / 2), v_0)
    print(f"Using the condition v_r(L/2) = v_0, we solve '{equation_for_k}' for k.")
    k_solution = sympy.solve(equation_for_k, k)[0]
    print(f"The constant k = {k_solution}")

    # Substitute k back into the general formula to get the specific velocity profile.
    v_r = v_r_general.subs(k, k_solution)
    print(f"The final river velocity profile is: v_r(y) = {v_r}\n")


    # Step 2: Calculate the downstream displacement for the outward journey.
    # The boat travels from y = 0 to y = 3L/4.
    # The time to travel a small distance dy is dt = dy / v.
    # The downstream drift during this time is dx = v_r(y) * dt = v_r(y) * (dy / v).
    print("Step 2: Calculate the downstream displacement for the outward journey (from y=0 to y=3L/4).")
    # To find the total downstream displacement, we integrate dx from y=0 to y=3L/4.
    x_outward = sympy.integrate(v_r / v, (y, 0, 3 * L / 4))
    print(f"The displacement is the integral of (v_r(y) / v) dy from 0 to 3L/4.")
    print(f"Displacement on outward journey = {x_outward}\n")

    # Step 3: Calculate the downstream displacement for the return journey.
    # The boat travels from y = 3L/4 back to y = 0.
    # The time taken to cover the same cross-stream distance dy is the same, so the drift is also the same.
    # The integral will be identical.
    print("Step 3: Calculate the downstream displacement for the return journey (from y=3L/4 to y=0).")
    print("The magnitude of time taken to cross any given band dy is the same, so the drift is the same.")
    x_return = x_outward
    print(f"Displacement on return journey = {x_return}\n")


    # Step 4: Sum the displacements to find the total displacement.
    print("Step 4: Sum the outward and return displacements for the total.")
    total_displacement = x_outward + x_return
    
    # Simplify and display the final result.
    final_expression = sympy.simplify(total_displacement)
    
    # The output from sympy is L*v_0*9/(8*v). Let's format it for clarity.
    num, den = final_expression.as_numer_denom()
    
    print("The final equation for the total downstream distance is:")
    # We manually format the print output to match the desired "each number" format.
    print("Distance = (9 * L * v_0) / (8 * v)")


solve_boat_problem()
<<< (9 * L * v_0) / (8 * v) >>>