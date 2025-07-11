import math

def solve_horse_water_problem(n, m):
    """
    Calculates the maximum amount of water left after a horse travels m km.

    Args:
        n (int): The number of 100-liter water loads at the start.
        m (float): The total distance to travel in km.
    """
    
    # The problem statement assumes n*100 > m, so the horse can always reach.
    # We can proceed with the calculation.

    initial_water = float(n * 100)
    distance_to_travel = float(m)
    
    water_left = initial_water
    
    # This list will hold the parts of our final equation.
    equation_parts = [f"{initial_water:.2f}"]
    
    # Iterate from the initial number of trips (n) down to 1.
    # Each iteration represents a segment of the journey with a constant consumption rate.
    for k in range(n, 0, -1):
        # If we've already covered the whole distance, stop.
        if distance_to_travel <= 1e-9:
            break

        # The consumption rate for this segment is (2k-1) liters per km.
        consumption_rate = 2 * k - 1
        
        # The maximum possible distance for a segment (where k > 1) is one where
        # exactly 100L are consumed, allowing us to reduce the number of trips.
        # If k=1, the horse just goes forward, so the segment can be infinitely long.
        max_dist_for_segment = 100.0 / consumption_rate if k > 1 else float('inf')
        
        # The actual distance we will travel in this segment is the smaller of
        # the remaining distance or the max possible segment distance.
        dist_in_this_segment = min(distance_to_travel, max_dist_for_segment)
        
        # Calculate the water consumed in this segment.
        water_consumed = consumption_rate * dist_in_this_segment
        
        # Update the total water left and the remaining distance.
        water_left -= water_consumed
        distance_to_travel -= dist_in_this_segment
        
        # Add this step to our equation string, showing each number.
        equation_parts.append(f"({consumption_rate} * {dist_in_this_segment:.2f})")

    # Format the final equation string for printing.
    final_equation = " - ".join(equation_parts)

    print("The maximum amount of water left is calculated as follows:")
    print(f"Final Water = {final_equation}")
    
    final_answer = water_left
    print(f"\nResult: The maximum amount of water left is {final_answer:.2f} liters.")
    
    # Return the final answer in the specified format.
    print(f"\n<<<{final_answer:.2f}>>>")

# --- User Input ---
# You can change these values to test different scenarios.
n_loads = 3  # Represents 3 * 100 = 300 liters of water
m_distance = 50 # Represents a 50 km journey

# Execute the function with the provided inputs
solve_horse_water_problem(n_loads, m_distance)