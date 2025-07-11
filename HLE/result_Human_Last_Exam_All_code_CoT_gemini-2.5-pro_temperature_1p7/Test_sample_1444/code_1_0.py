import math

def solve_tree_journey():
    """
    Calculates the total distance and direction of the walking tree's journey.
    """
    # Given data
    initial_position_degrees = 45.0
    position_after_first_walk_degrees = 44.8
    final_position_degrees = 34.4
    first_walk_distance_meters = 100

    # Determine the direction of travel
    # The angle to Polaris represents latitude. A decrease means moving south.
    direction = "South"

    # --- Calculations ---

    # 1. Calculate the angular change during the first walk
    first_walk_angular_change = initial_position_degrees - position_after_first_walk_degrees

    # 2. Calculate the total angular change over the entire period
    total_angular_change = initial_position_degrees - final_position_degrees

    # 3. Calculate the total distance traveled
    # We use the first walk to calibrate the distance per degree of change.
    # Total Distance = (Total Angular Change / First Walk's Angular Change) * First Walk's Distance
    total_distance_meters = (total_angular_change / first_walk_angular_change) * first_walk_distance_meters
    total_distance_km = total_distance_meters / 1000

    # --- Output Results ---
    
    print("Direction of travel:")
    print(f"The tree's viewing angle to Polaris decreased from {initial_position_degrees}° to {final_position_degrees}°, meaning its latitude decreased. Therefore, the tree walked {direction}.")
    
    print("\nTotal distance calculation:")
    # Per the instructions, showing each number in the final equation
    print(f"The total distance is calculated by scaling the total angular change ({total_angular_change:.1f}°) with the reference from the first walk ({first_walk_angular_change:.1f}° for {first_walk_distance_meters}m).")
    print(f"Equation: Total Distance = ({total_angular_change:.1f} / {first_walk_angular_change:.1f}) * {first_walk_distance_meters}")
    
    print(f"\nApproximate total distance traveled: {total_distance_km:.1f} km")

    # --- Final Answer Formatting ---
    
    # Final answer as per instruction: Nearest Integer(Total Distance in km * 10)
    final_answer = round(total_distance_km * 10)
    # The final answer is wrapped as requested
    print(f"\n<<<>>>")

if __name__ == '__main__':
    solve_tree_journey()