import sympy

def solve_boat_problem():
    """
    Calculates the total downstream drift of a boat crossing a river,
    turning around at 3/4 of the way, and returning to the starting bank.
    """
    # Step 1: Define all symbolic variables for the problem.
    L, v_0, v, y, k = sympy.symbols('L v_0 v y k', real=True, positive=True)
    print("--- Step 1: Define physical quantities as symbolic variables ---")
    print(f"River width: {L}")
    print(f"Maximum flow velocity at the center: {v_0}")
    print(f"Boat speed relative to water (perpendicular to flow): {v}")
    print(f"Variable for distance from the starting bank: {y}\n")

    # Step 2: Determine the river's flow velocity function, v_flow(y).
    # The flow profile is parabolic: v_flow = k*y*(L-y).
    # We find 'k' by using the condition that velocity is v_0 at y = L/2.
    # Equation: v_0 = k * (L/2) * (L - L/2)
    k_val = sympy.solve(k * (L/2) * (L - L/2) - v_0, k)[0]
    v_flow = k_val * y * (L - y)
    print("--- Step 2: Determine the river's flow velocity profile ---")
    print(f"The velocity profile is parabolic: v_flow(y) = k*y*(L-y)")
    print(f"Using the condition v_flow(L/2) = v_0, we find k = {k_val}")
    print(f"The complete flow velocity function is: v_flow(y) = {v_flow}\n")

    # Step 3: Calculate the downstream drift for the outward journey.
    # The boat travels from y=0 to y=3L/4.
    # The boat's velocity across the river is 'v'.
    # The drift (x_out) is the integral of (v_flow / v) with respect to y.
    # The limits of integration are from 0 to 3L/4.
    print("--- Step 3: Calculate the drift during the outward journey (from y=0 to y=3L/4) ---")
    x_out = sympy.integrate(v_flow / v, (y, 0, 3*L/4))
    print(f"The outward drift, x_out, is the integral of (v_flow(y) / v) dy from 0 to 3L/4.")
    print(f"x_out = {sympy.simplify(x_out)}\n")

    # Step 4: Calculate the downstream drift for the return journey.
    # The boat travels from y=3L/4 back to y=0.
    # The integral for the return drift, x_back, is mathematically identical
    # to the one for the outward drift.
    print("--- Step 4: Calculate the drift during the return journey (from y=3L/4 to y=0) ---")
    x_back = sympy.integrate(v_flow / v, (y, 0, 3*L/4))
    print(f"The return drift, x_back, is calculated similarly and is equal to the outward drift.")
    print(f"x_back = {sympy.simplify(x_back)}\n")

    # Step 5: Sum the drifts to find the total distance.
    total_drift = x_out + x_back
    final_simplified_drift = sympy.simplify(total_drift)
    print("--- Step 5: Calculate the total downstream distance ---")
    print("The total distance is the sum of the outward and return drifts.")
    print(f"Total Distance = x_out + x_back = {sympy.simplify(x_out)} + {sympy.simplify(x_back)}")

    # Final Answer Output
    # The numbers in the final equation are 9 and 8.
    print("\n--- Final Answer ---")
    print("The total distance between the boat's returning position and its original starting point is:")
    print(f"Total Distance = {final_simplified_drift}")

if __name__ == '__main__':
    solve_boat_problem()