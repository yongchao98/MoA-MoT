import math

def solve_horse_problem(n, m):
    """
    Calculates the maximum amount of water left after a horse travels m km.

    Args:
        n (int): The number of 100-liter water containers at the origin.
        m (int): The distance to the destination in km.
    """
    # Validate input: Ensure there's theoretically enough water to reach the destination.
    if n * 100.0 <= m:
        print("Error: The total initial water is not enough to travel the distance, even with a perfect 1 L/km consumption rate.")
        print("<<<0.0000>>>")
        return

    # Keep track of the distance covered and build the summation part of the equation
    cumulative_dist = 0.0
    sum_dist_str_list = []

    # Loop through the stages. Each stage i corresponds to depleting one container of water.
    for i in range(1, n):
        # Number of containers being moved in the current stage
        num_containers = n - (i - 1)
        
        # Consumption rate is (2k-1) where k is the number of containers
        consumption_rate = 2 * num_containers - 1
        
        # The maximum distance covered in a stage by consuming one 100L container
        dist_for_one_container = 100.0 / consumption_rate

        # Check if the destination m falls within the current stage
        if cumulative_dist + dist_for_one_container >= m:
            
            # The amount of water available at the start of this stage
            water_at_stage_start = num_containers * 100.0
            
            # Distance that needs to be traveled in this final leg
            dist_to_travel_in_this_stage = m - cumulative_dist
            
            # Water consumed in this final leg
            water_consumed = dist_to_travel_in_this_stage * consumption_rate
            
            # Final amount of water left
            final_water = water_at_stage_start - water_consumed
            
            # --- Construct the output strings ---
            # Create the summation string for previous stages. It's '0' if this is the first stage.
            sum_str = " + ".join(sum_dist_str_list) if sum_dist_str_list else "0"
            
            print("The destination is reached in a stage where multiple trips are needed to move the water forward.")
            print(f"In this stage, the horse is moving {num_containers} containers, leading to a consumption rate of {consumption_rate} L/km.")
            
            print("\nThe final amount of water is the water at the start of this stage minus the water consumed to travel the remaining distance.")
            print("\nFinal Equation (water left = water_at_stage_start - (distance_in_stage * rate)):")
            print(f"Maximum water left = {num_containers}*100 - ({m} - ({sum_str})) * {consumption_rate}")
            
            print(f"\nResult: {final_water:.4f} liters")
            print(f"\n<<<{final_water:.4f}>>>")
            return
            
        # If the destination is not in this stage, complete the stage and update trackers.
        cumulative_dist += dist_for_one_container
        # The string term uses num_containers, which corresponds to the rate for this stage
        sum_dist_str_list.append(f"100/(2*{num_containers}-1)")

    # If the loop finishes, the horse is on its last 100L container.
    # The water left at the start of this final stage is exactly 100L.
    water_at_stage_start = 100.0
    consumption_rate = 1.0  # Only one-way travel now
    
    # Distance to travel in this final one-way trip
    dist_to_travel_in_this_stage = m - cumulative_dist
    
    # Water consumed in this one-way trip
    water_consumed = dist_to_travel_in_this_stage * consumption_rate
    
    # Final amount of water left
    final_water = water_at_stage_start - water_consumed

    # --- Construct the output strings ---
    sum_str = " + ".join(sum_dist_str_list)
    
    print("The destination is far enough that the horse first creates depots until only 100L of water remains.")
    print("Then, it travels the final leg of the journey with that last container.")

    print("\nThe final amount of water is 100L minus the water consumed during the final one-way trip.")
    print("\nFinal Equation (water left = 100 - (distance_in_stage * 1)):")
    # The term (m - sum(...)) is the distance_in_stage
    print(f"Maximum water left = 100 - ({m} - ({sum_str}))")
    
    print(f"\nResult: {final_water:.4f} liters")
    print(f"\n<<<{final_water:.4f}>>>")

# --- Example Usage ---
# You can change the values of n and m to solve for different scenarios.
n_containers = 20 # example: 20 * 100 = 2000 liters
distance_km = 1500 # example: 1500 km

# n_containers = 3
# distance_km = 70

solve_horse_problem(n=n_containers, m=distance_km)