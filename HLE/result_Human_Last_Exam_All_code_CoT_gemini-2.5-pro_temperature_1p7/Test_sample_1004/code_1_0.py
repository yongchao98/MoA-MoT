from fractions import Fraction

def solve_boat_problem():
    """
    This function calculates the final downstream displacement of the boat.
    The result is symbolic, expressed as a coefficient multiplied by (v_0 * L / v).
    """

    # Step 1: Define coefficients for the drift calculation based on the physics.
    # The journey is broken into two parts: out and back. The drift is identical for both parts.
    # Total Drift D = 2 * (Drift for one way)
    # Drift for one way = Integral of v_flow(y)/v dy from y=0 to y=3L/4.
    # We split this integral at y=L/2.

    # Part 1 of the integral (from y=0 to y=L/2):
    # Integral results in a drift of (1/4) * (v_0 * L / v)
    coeff1 = Fraction(1, 4)

    # Part 2 of the integral (from y=L/2 to y=3L/4):
    # Integral results in a drift of (3/16) * (v_0 * L / v)
    coeff2 = Fraction(3, 16)

    # Step 2: Calculate the total coefficient for one-way drift.
    one_way_drift_coeff = coeff1 + coeff2

    # Step 3: Calculate the total coefficient for the round trip.
    total_drift_coeff = 2 * one_way_drift_coeff
    
    # Step 4: Print the derivation and the final formula.
    print("Let D be the distance between the returning position and the starting point.")
    print("The final distance D is a product of a numerical coefficient and the term (v_0 * L / v).\n")
    print("The derivation is as follows:")
    
    # Explain the formula structure
    print("D = 2 * (Drift_out) = 2 * (Drift_0_to_L/2 + Drift_L/2_to_3L/4)")
    print("The calculation for the coefficient of (v_0 * L / v) is:")
    
    # Show the numerical calculation using the coefficients
    print(f"D_coeff = 2 * ({coeff1} + {coeff2})")
    print(f"D_coeff = 2 * ({one_way_drift_coeff})")
    print(f"D_coeff = {total_drift_coeff}")
    
    # Display the final formula
    num = total_drift_coeff.numerator
    den = total_drift_coeff.denominator
    
    print("\nThus, the final distance D is given by the formula:")
    print(f"D = ({num}/{den}) * (v_0 * L / v)")

# Run the solver
solve_boat_problem()