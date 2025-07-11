import fractions

def calculate_boat_displacement():
    """
    This function calculates the final downstream displacement of a boat
    crossing a river with a parabolic flow profile.

    The final distance is expressed as a coefficient multiplied by (v_0 * L / v).
    This script calculates that coefficient.
    """

    # The calculation for the downstream drift on the outbound trip
    # (from y=0 to y=3L/4) yields a coefficient of (9/16).
    # x_out = (9/16) * (v_0 * L / v)
    coeff_out = fractions.Fraction(9, 16)

    # The calculation for the downstream drift on the return trip
    # (from y=3L/4 to y=0) yields the same coefficient.
    # x_back = (9/16) * (v_0 * L / v)
    coeff_back = fractions.Fraction(9, 16)

    # The total displacement is the sum of the two drifts.
    total_coeff = coeff_out + coeff_back

    # Extract the numerator and denominator for the final output.
    numerator = total_coeff.numerator
    denominator = total_coeff.denominator

    # Print the final equation showing each number.
    print("The final equation for the total distance is:")
    print(f"Distance = ({numerator}/{denominator}) * (v_0 * L / v)")

# Run the calculation and print the result.
calculate_boat_displacement()