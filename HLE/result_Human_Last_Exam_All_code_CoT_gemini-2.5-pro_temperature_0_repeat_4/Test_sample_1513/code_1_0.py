def find_remote_caldera_town():
    """
    This function identifies an island town based on specific geographical and naming criteria.
    """
    # Step 1: Define the criteria from the riddle.
    min_distance_from_inhabited_land_km = 500
    geological_feature = "volcanic caldera bay"
    shared_name_feature = "town and bay"

    # Step 2: Identify the location that matches the criteria from a knowledge base.
    town_name = "Hanga Roa"
    bay_name = "Hanga Roa"
    island_name = "Easter Island (Rapa Nui)"
    actual_distance_to_nearest_inhabited_km = 2075 # Distance to Pitcairn Islands

    # Step 3: Verify that the location meets all conditions.
    if (town_name == bay_name and
        actual_distance_to_nearest_inhabited_km > min_distance_from_inhabited_land_km):

        # Step 4: Print the results, including the number from the prompt.
        print(f"The island town is {town_name}.")
        print(f"It is located on {island_name}.")
        print(f"The town sits on {bay_name} Bay, which is part of a bay formed by volcanic activity.")
        print(f"The town and the bay share the name '{town_name}'.")
        print(f"The island is more than {min_distance_from_inhabited_land_km} kilometers from another inhabited island.")
        print(f"Its actual distance to the nearest inhabited land (Pitcairn Islands) is approximately {actual_distance_to_nearest_inhabited_km} km.")

find_remote_caldera_town()
<<<Hanga Roa>>>