def solve_river_problem():
    """
    This script calculates the total downstream drift of a boat crossing a river.
    It provides a step-by-step derivation of the final formula.
    """

    print("Step-by-step derivation of the total drift:")
    print("-" * 50)

    # Step 1: Define the river velocity profile v_x(y)
    print("1. The river's flow velocity v_x as a function of distance y from the shore is a triangular profile:")
    print("   - v_x(y) = (2 * v_0 / L) * y, for 0 <= y <= L/2")
    print("   - v_x(y) = (2 * v_0 / L) * (L - y), for L/2 < y <= L")
    print("\n")

    # Step 2: Formulate the drift integral
    print("2. The drift 'dx' is v_x * dt. The boat's cross-stream speed is v_y = dy/dt = v.")
    print("   Therefore, dt = dy / v.")
    print("   The drift can be found by integrating: Drift = integral(v_x(y) / v) dy.")
    print("\n")

    # Step 3: Calculate the outbound drift (from y=0 to y=3L/4)
    print("3. Calculate the outbound drift, which must be split into two integrals at y = L/2:")
    print("   Drift_out = (1/v) * [ integral from 0 to L/2 of (2*v_0/L)*y dy + integral from L/2 to 3L/4 of (2*v_0/L)*(L-y) dy ]")
    print("   - First integral (0 to L/2): (2*v_0 / (L*v)) * [y^2/2]_0^(L/2) = (v_0*L) / (4*v)")
    print("   - Second integral (L/2 to 3L/4): (2*v_0 / (L*v)) * [L*y - y^2/2]_(L/2)^(3L/4) = (3*v_0*L) / (16*v)")
    print("   Total Outbound Drift = (v_0*L)/(4*v) + (3*v_0*L)/(16*v) = (7 * v_0 * L) / (16 * v)")
    print("\n")

    # Step 4: Calculate the return drift (from y=3L/4 to y=0)
    print("4. The return trip covers the same path in reverse.")
    print("   The integral for the return drift is identical to the outbound drift.")
    print("   Total Return Drift = (7 * v_0 * L) / (16 * v)")
    print("\n")

    # Step 5: Calculate the total drift
    print("5. The total drift is the sum of the outbound and return drifts.")
    print("   Total Drift = Drift_out + Drift_return")
    print("   Total Drift = (7*v_0*L)/(16*v) + (7*v_0*L)/(16*v)")
    print("   Total Drift = (14 * v_0 * L) / (16 * v)")
    print("\n")

    # Step 6: Final Answer
    print("6. Simplifying the fraction gives the final answer:")
    numerator = 7
    denominator = 8
    print(f"Total Drift = ({numerator} * v_0 * L) / ({denominator} * v)")

solve_river_problem()