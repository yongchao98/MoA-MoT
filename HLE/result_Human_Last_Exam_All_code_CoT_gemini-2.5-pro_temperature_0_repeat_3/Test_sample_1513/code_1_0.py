def solve_geography_riddle():
    """
    This script solves a geographical riddle by encoding the clues
    and presenting the location that matches them.
    """

    # --- Clues from the riddle ---
    # The final equation is the combination of these facts leading to the answer.
    # Fact 1: Distance requirement
    min_distance_km = 500

    # Fact 2: Geological features
    geological_feature = "volcanic caldera"
    location_type = "bay"

    # Fact 3: Naming convention
    # The town's name is the same as the bay's name.

    # --- The Solution ---
    town = "Hanga Roa"
    island = "Easter Island (Rapa Nui)"
    bay = "Hanga Roa Bay"
    distance_to_nearest_inhabited_land_km = 2075 # Distance to Pitcairn Islands

    # --- Presenting the Answer ---
    print(f"Solving the riddle with the following facts:")
    print(f"1. The town must be more than {min_distance_km} km from another inhabited island.")
    print(f"2. The town sits on a {location_type} on a volcanic island with a {geological_feature}.")
    print(f"3. The town and the bay share the same name.")
    print("-" * 20)

    print(f"The island town is: {town}")
    print(f"It is the capital of {island}, a volcanic island.")
    print(f"It sits on {bay}, fulfilling the shared name requirement.")
    print(f"The nearest inhabited land is over {distance_to_nearest_inhabited_land_km} km away, which is greater than the required {min_distance_km} km.")

solve_geography_riddle()