import fractions

def solve_river_crossing_problem():
    """
    Calculates the final drift distance for the boat and prints the derivation.
    """
    print("Solving the river boat problem step-by-step.")
    print("Let L be the river width, v0 the max flow speed, and v the boat's speed relative to the water.")
    print("="*60)

    print("\nStep 1: Define the River Flow Velocity Profile")
    print("The flow velocity v_flow(y) at a distance y from the home bank is modeled as a")
    print("piecewise linear function:")
    print("  - For 0 <= y <= L/2:   v_flow(y) = (2 * v0 / L) * y")
    print("  - For L/2 <= y <= L:   v_flow(y) = (2 * v0 / L) * (L - y)")
    print("This profile has v_flow=0 at the banks (y=0, y=L) and v_flow=v0 at the center (y=L/2).")

    print("\nStep 2: Calculate Downstream Drift for the Outbound Journey (y=0 to y=3L/4)")
    print("The boat's travel time is related to its y-position by dt = dy/v.")
    print("The downstream drift is the integral of v_flow(y) * dt.")
    print("Drift = integral(v_flow(y) / v dy) from y=0 to y=3L/4.")
    print("We split the integral at the river's center, y=L/2.")

    # The drift for any part of the journey is proportional to (v0*L/v).
    # We will calculate the dimensionless coefficient.

    # Part A: Drift from y=0 to y=L/2
    # The integral of ( (2*v0/L)*y / v ) dy from 0 to L/2 gives (1/4)*(v0*L/v)
    coeff_A = fractions.Fraction(1, 4)
    print(f"\n  - Part A (y=0 to L/2): The drift is ({coeff_A.numerator}/{coeff_A.denominator}) * (v0*L/v).")

    # Part B: Drift from y=L/2 to y=3L/4
    # The integral of ( (2*v0/L)*(L-y) / v ) dy from L/2 to 3L/4 gives (3/16)*(v0*L/v)
    coeff_B = fractions.Fraction(3, 16)
    print(f"  - Part B (y=L/2 to 3L/4): The drift is ({coeff_B.numerator}/{coeff_B.denominator}) * (v0*L/v).")

    # Total Outbound Drift
    total_outbound_coeff = coeff_A + coeff_B
    print(f"\n  The total outbound drift is the sum of these parts:")
    print(f"  ({coeff_A}) * (v0*L/v) + ({coeff_B}) * (v0*L/v) = ({total_outbound_coeff}) * (v0*L/v).")

    print("\nStep 3: Calculate Downstream Drift for the Return Journey (y=3L/4 to y=0)")
    print("The boat travels back across the same current profile for the same amount of time.")
    print("Therefore, the drift during the return journey is identical to the outbound journey.")
    return_coeff = total_outbound_coeff
    print(f"  The return drift is also ({return_coeff}) * (v0*L/v).")

    print("\nStep 4: Calculate the Total Distance from the Starting Point")
    print("The total distance is the sum of the outbound and return drifts.")
    total_coeff = total_outbound_coeff + return_coeff
    print(f"  Total Distance = ({total_outbound_coeff}) * (v0*L/v) + ({return_coeff}) * (v0*L/v)")
    print(f"                   = ({total_coeff}) * (v0*L/v)")


    print("\n" + "="*60)
    print("Final Answer")
    print("="*60)
    print("The final distance between the boat's starting point and its returning position is:")
    print(f"\n    Distance = ({total_coeff.numerator} / {total_coeff.denominator}) * (v0 * L / v)\n")

    print("The numbers in the final equation's fraction are:")
    print(f"  Numerator:   {total_coeff.numerator}")
    print(f"  Denominator: {total_coeff.denominator}")

if __name__ == "__main__":
    solve_river_crossing_problem()