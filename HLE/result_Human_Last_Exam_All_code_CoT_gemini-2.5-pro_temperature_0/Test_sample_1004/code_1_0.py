def solve_river_crossing_problem():
    """
    This function provides a step-by-step derivation for the river crossing problem
    and prints the final answer.
    """

    print("This script solves for the total downstream drift of the boat.")
    print("Let L be the river width, v_0 be the max flow velocity, and v be the boat's speed relative to the water.")
    print("-" * 60)

    # Step 1: Derivation of the river's velocity profile u(y)
    print("Step 1: Determine the river's velocity profile u(y).")
    print("The flow velocity u(y) is proportional to the distance from both shores (y and L-y).")
    print("This can be modeled as: u(y) = k * y * (L - y).")
    print("Given u(L/2) = v_0, we find k = 4 * v_0 / L^2.")
    print("So, the velocity profile is: u(y) = (4 * v_0 / L^2) * (L*y - y^2).")
    print("-" * 60)

    # Step 2: Calculation of the outbound drift (x1)
    print("Step 2: Calculate the downstream drift for the outbound trip (from y=0 to y=3L/4).")
    print("The boat's velocity across the river is dy/dt = v, so dt = dy/v.")
    print("The drift is the integral of the river's velocity over time: x1 = integral(u(y) dt).")
    print("Substituting dt, we get: x1 = integral from 0 to 3L/4 of (u(y)/v) dy.")
    print("The integral of (4*v_0/(v*L^2)) * (L*y - y^2) is (4*v_0/(v*L^2)) * [L*y^2/2 - y^3/3].")
    print("Evaluating this from 0 to 3L/4 yields: x1 = (9 * v_0 * L) / (16 * v).")
    print("-" * 60)

    # Step 3: Calculation of the inbound drift (x2)
    print("Step 3: Calculate the downstream drift for the inbound trip (from y=3L/4 to y=0).")
    print("The boat's velocity across the river is dy/dt = -v, so dt = dy/(-v).")
    print("The drift is: x2 = integral from 3L/4 to 0 of (u(y)/(-v)) dy.")
    print("This integral evaluates to the same result as the outbound trip: x2 = (9 * v_0 * L) / (16 * v).")
    print("-" * 60)

    # Step 4: Calculation of the total drift
    print("Step 4: Calculate the total downstream drift.")
    print("Total Drift = x1 + x2")
    print("Total Drift = (9 * v_0 * L) / (16 * v) + (9 * v_0 * L) / (16 * v)")
    print("Total Drift = 2 * (9 * v_0 * L) / (16 * v)")
    print("-" * 60)

    # Final Answer
    print("The final simplified distance between the starting and returning points is:")
    
    numerator_coeff = 9
    denominator_coeff = 8
    
    # Outputting the final equation with its numerical components
    print(f"({numerator_coeff} * v_0 * L) / ({denominator_coeff} * v)")

# Run the solver
solve_river_crossing_problem()