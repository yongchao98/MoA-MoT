def calculate_boat_drift():
    """
    Calculates the total downstream drift of a boat based on a derived formula.

    The boat crosses a river of width L, where the flow velocity is parabolic,
    reaching a maximum of v0 at the center. The boat's speed relative to the
    water is v. It travels 3/4 of the way across, then immediately returns
    to the starting bank.

    The final derived formula for the total drift distance is:
    Total Drift = (9 * v0 * L) / (8 * v)
    """

    # --- Define example problem parameters ---
    # River width in meters
    L = 100.0
    # Maximum flow velocity at the river's center in m/s
    v0 = 2.0
    # Boat's speed relative to the water (perpendicular to flow) in m/s
    v = 5.0

    print("--- River Crossing Problem ---")
    print("This script calculates the final displacement of the boat.")
    print("\nGiven Parameters:")
    print(f"  River width (L): {L} m")
    print(f"  Max flow velocity (v0): {v0} m/s")
    print(f"  Boat's relative speed (v): {v} m/s")

    # --- Calculate the final drift distance using the formula ---
    
    # Numerator of the formula: 9 * v0 * L
    numerator = 9 * v0 * L
    
    # Denominator of the formula: 8 * v
    denominator = 8 * v
    
    total_drift = numerator / denominator

    print("\nCalculation:")
    print("The formula for the total downstream drift is: (9 * v0 * L) / (8 * v)")
    print("Substituting the given values into the equation:")
    # Here we output each number in the final equation as requested
    print(f"  Total Drift = (9 * {v0} * {L}) / (8 * {v})")
    print(f"  Total Drift = {numerator} / {denominator}")
    
    print("\nResult:")
    print(f"The distance between the boat's returning position and its original starting point is {total_drift:.2f} meters.")

# Execute the function
calculate_boat_drift()