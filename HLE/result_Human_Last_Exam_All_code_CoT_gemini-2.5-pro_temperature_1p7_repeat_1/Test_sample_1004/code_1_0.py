def calculate_boat_displacement():
    """
    Calculates the total downstream displacement of a boat crossing a river and returning.

    The problem parameters are:
    L: width of the river
    v0: maximum flow velocity at the center of the river
    v: the boat's speed relative to the water, perpendicular to the flow

    The final formula for the total downstream displacement is (7 * v0 * L) / (8 * v).
    """
    # Assigning example values to the variables for demonstration
    L = 100.0  # River width in meters
    v0 = 5.0   # Maximum flow velocity at the center in m/s
    v = 4.0    # Boat's speed relative to the water in m/s

    print("This script calculates the final downstream displacement of the boat.")
    print("Given values for demonstration:")
    print(f"L (River width) = {L} m")
    print(f"v0 (Max flow velocity) = {v0} m/s")
    print(f"v (Boat's relative speed) = {v} m/s")
    print("-" * 40)

    # The total downstream displacement is given by the derived formula.
    # The script will now use the formula to calculate the result and
    # display the intermediate steps of the calculation.

    numerator = 7 * v0 * L
    denominator = 8 * v
    total_displacement = numerator / denominator

    print("The derived formula for the total downstream displacement is:")
    print("Displacement = (7 * v0 * L) / (8 * v)")
    print("\nSubstituting the given values into the formula:")

    # Outputting each number in the final equation as requested
    print(f"Displacement = (7 * {v0} * {L}) / (8 * {v})")
    print(f"             = {numerator} / {denominator}")
    print(f"             = {total_displacement:.4f} meters")
    print("-" * 40)
    print("The distance between the boat's returning position on the home bank and its original starting point is {:.4f} meters.".format(total_displacement))


calculate_boat_displacement()
