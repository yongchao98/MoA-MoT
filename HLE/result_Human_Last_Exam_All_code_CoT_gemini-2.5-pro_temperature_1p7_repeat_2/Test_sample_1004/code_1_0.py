def calculate_boat_displacement():
    """
    Calculates the downstream displacement of a boat based on a derived formula.

    The problem describes a boat crossing a river with a parabolic flow profile,
    traveling 3/4 of the way across, and then returning to the starting bank.

    The derivation involves integrating the river's velocity against the time
    the boat spends at each cross-stream position 'y'.

    The river velocity is v_flow(y) = (4*v0/L^2) * y * (L-y).
    The boat's cross-stream velocity is v.

    The total downstream displacement is the sum of the displacement on the
    outward and return journeys.
    - Drift out (y=0 to 3L/4): integral of v_flow(y) * dt, where dt = dy/v.
    - Drift back (y=3L/4 to 0): integral of v_flow(y) * dt, where dt = dy/(-v).

    Both integrals evaluate to the same value: (9 * v0 * L) / (16 * v).
    The total displacement is the sum of these two, which simplifies to the
    final formula used in this function.
    """

    # The final derived formula for the total downstream displacement is:
    # Displacement = (9 / 8) * (L * v0) / v

    print("The final equation for the displacement is derived by integrating the boat's downstream drift over its entire journey.")
    print("Final Formula: Displacement = (numerator / denominator) * (L * v0) / v")
    
    # As requested, here are the numbers in the final equation:
    numerator = 9
    denominator = 8
    
    print(f"\nThe number in the numerator is: {numerator}")
    print(f"The number in the denominator is: {denominator}")
    
    print(f"\nSo, the final equation is: Displacement = ({numerator} / {denominator}) * (L * v0) / v")
    print("\nLet's calculate the displacement with some example values:")

    # Example values for L, v0, and v
    L = 100.0  # meters
    v0 = 2.0   # m/s (max flow velocity at center)
    v = 4.0    # m/s (boat's speed relative to water)

    print(f"  River width L = {L} m")
    print(f"  Max flow velocity v0 = {v0} m/s")
    print(f"  Boat's relative speed v = {v} m/s")

    # Calculate the displacement using the formula
    displacement = (numerator / denominator) * (L * v0) / v

    print(f"\nThe calculated distance between the start and end points is: {displacement} meters")

if __name__ == '__main__':
    calculate_boat_displacement()