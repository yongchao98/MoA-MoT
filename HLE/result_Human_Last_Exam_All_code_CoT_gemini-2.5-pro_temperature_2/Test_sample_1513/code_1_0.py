def solve_riddle():
    """
    This script provides the answer to the geographical riddle by breaking down the requirements.
    """
    # Define the elements of the answer based on the riddle's clues.
    town_name = "Hanga Roa"
    bay_name = "Hanga Roa Bay"
    geological_feature = "volcanic caldera"
    required_distance_km = 500
    actual_distance_to_neighbor_km = 2075 # Approximate distance to Pitcairn Islands

    # Print the solution, explaining how each condition is met.
    print(f"The island town that matches the description is: {town_name}")
    print(f"1. Name: The town '{town_name}' sits on '{bay_name}', so the names are shared.")
    print(f"2. Geology: The town is on a bay formed by the island's volcanic landscape, which includes a major '{geological_feature}'.")
    print(f"3. Remoteness: The distance to the nearest inhabited island is ~{actual_distance_to_neighbor_km} km, which is more than the required {required_distance_km} km.")

solve_riddle()