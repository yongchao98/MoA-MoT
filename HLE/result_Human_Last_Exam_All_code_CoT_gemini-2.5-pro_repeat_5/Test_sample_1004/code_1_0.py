import sys

def solve_boat_problem(L, v0, v):
    """
    Calculates the final downstream displacement of the boat.

    Args:
        L (float): Width of the river.
        v0 (float): Maximum flow velocity of the river at the center.
        v (float): The boat's speed relative to the water, perpendicular to the flow.

    Returns:
        float: The total downstream distance from the starting point.
    """
    if v == 0:
        print("Error: The boat's relative speed 'v' cannot be zero.", file=sys.stderr)
        return None
    
    # The final derived formula for the total downstream distance is:
    # distance = (7 * v0 * L) / (8 * v)
    distance = (7 * v0 * L) / (8 * v)
    
    # We print the final equation with the values substituted in
    print("Problem Parameters:")
    print(f"River width (L) = {L} m")
    print(f"Max river velocity (v0) = {v0} m/s")
    print(f"Boat's relative speed (v) = {v} m/s")
    print("-" * 30)
    print("Final Calculation:")
    print(f"Distance = (7 * v0 * L) / (8 * v)")
    print(f"Distance = (7 * {v0} * {L}) / (8 * {v})")
    print(f"Distance = {distance:.2f} m")
    
    return distance

# --- Example Usage ---
# You can change these values to see the result for different scenarios.
river_width = 100.0  # meters
max_flow_velocity = 3.0   # m/s
boat_relative_speed = 5.0 # m/s

# Calculate and print the result
solve_boat_problem(river_width, max_flow_velocity, boat_relative_speed)