import math

def solve_walking_tree():
    """
    Calculates the total distance traveled and direction of the Ancient Walking Tree.
    """

    # 1. Define the known values from the problem description.
    initial_latitude = 45.0  # in degrees
    final_latitude = 34.4  # in degrees
    
    # A standard approximation for the length of one degree of latitude on Earth.
    # The problem asks for an approximate distance, so this is a reasonable constant to use.
    km_per_degree_latitude = 111.2

    # 2. Determine the direction.
    direction = "South"
    print(f"The tree's latitude decreased, so it walked {direction}.")

    # 3. Calculate the total distance traveled.
    # The total distance is the total change in latitude multiplied by the conversion factor.
    total_latitude_change = initial_latitude - final_latitude
    total_distance_km = total_latitude_change * km_per_degree_latitude

    print("\nCalculating the total distance:")
    print("The equation for total distance is: (Initial Latitude - Final Latitude) * Kilometers per Degree")
    print(f"Total Distance (km) = ({initial_latitude} - {final_latitude}) * {km_per_degree_latitude}")
    print(f"Total Distance (km) = {total_latitude_change:.1f} * {km_per_degree_latitude}")
    print(f"Total Distance (km) = {total_distance_km:.2f}")

    # 4. Calculate the final answer in the requested format.
    final_answer = round(total_distance_km * 10)
    
    print("\nFinal requested value calculation:")
    print("Nearest Integer(Total Distance in km * 10)")
    print(f"Nearest Integer({total_distance_km:.2f} * 10)")
    print(f"Nearest Integer({total_distance_km * 10:.2f})")
    print(f"Final Answer Value: {final_answer}")
    
    return final_answer

# Run the calculation and store the final answer.
final_answer_value = solve_walking_tree()
print(f"<<<{final_answer_value}>>>")