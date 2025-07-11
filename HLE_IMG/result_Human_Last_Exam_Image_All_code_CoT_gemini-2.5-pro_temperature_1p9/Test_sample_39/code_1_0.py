def solve_bee_foraging():
    """
    Determines the bees' foraging location based on the waggle dance
    observed in the provided image.
    """

    # Step 1: The problem states that directions on the honeycomb frame
    # correspond to specific foraging locations. Let's define this map.
    foraging_map = {
        "North (up)": "In the Ambrosia apple orchard.",
        "West (left)": "In the McIntosh apple orchard.",
        "South (down)": "In the strawberry field.",
        "East (right)": "In the squash field."
    }

    print("Step 1: Establish the map of foraging locations.")
    print("-------------------------------------------------")
    for direction, location in foraging_map.items():
        print(f"Direction {direction} -> {location}")
    print("\n")

    # Step 2: Observe the image to find the waggle dance.
    # In the center-right of the frame, one bee is blurred from rapid
    # side-to-side movement. This is the "waggle dance." Surrounding
    # bees are oriented towards it, paying attention. The bee's dance
    # is oriented horizontally, pointing to the right.
    observed_direction = "East (right)"
    
    print("Step 2: Analyze the waggle dance in the image.")
    print("----------------------------------------------")
    print(f"A bee is performing a waggle dance, indicating a food source.")
    print(f"The direction of the dance on the comb is horizontal, pointing to the right.")
    print(f"This corresponds to the direction: {observed_direction}.\n")

    # Step 3: Match the observed direction to the foraging map to find the answer.
    final_location = foraging_map[observed_direction]

    print("Step 3: Conclude the most likely foraging location.")
    print("-----------------------------------------------------")
    print(f"Based on the dance direction, the bees are most likely foraging at the following location:")
    print(f"--> {final_location}")


solve_bee_foraging()