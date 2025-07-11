import sys

def solve_bee_foraging_location():
    """
    This script determines the bees' foraging location based on the waggle dance.
    """
    # Step 1: Define the known locations and directional angles.
    # North is 0 degrees, East is 90, South is 180, West is 270.
    print("Step 1: Establishing known information")
    crop_locations = {
        "North": "Ambrosia apple orchard",
        "West": "McIntosh apple orchard",
        "South": "strawberry field",
        "East": "large squash field"
    }
    cardinal_angles = {
        "North": 0,
        "East": 90,
        "South": 180,
        "West": 270
    }
    print("Crop locations and their directions from the hive are defined.\n")

    # Step 2: Analyze the visual and contextual clues from the problem.
    print("Step 2: Analyzing the clues")
    # In the image, the bee's dance is horizontal to the right, which is a 90-degree clockwise
    # turn from the vertical 'up' direction on the comb.
    dance_angle_from_vertical = 90
    print(f"The waggle dance in the image is at an angle of {dance_angle_from_vertical} degrees clockwise from vertical.")

    # The check happens 'before breakfast', implying it's morning. The sun is in the East.
    sun_direction_name = "East"
    sun_direction_angle = cardinal_angles[sun_direction_name]
    print(f"The time is morning, so the sun is in the {sun_direction_name} ({sun_direction_angle} degrees).\n")

    # Step 3: Apply the science of the waggle dance to calculate the food source direction.
    print("Step 3: Calculating the foraging direction")
    # On a vertical comb, the 'up' direction represents the sun's direction.
    # The food source direction is the sun's direction plus the angle of the dance.
    food_source_angle = (sun_direction_angle + dance_angle_from_vertical) % 360
    
    print("The final direction is calculated with the equation:")
    print(f"Sun's Angle + Dance's Angle = Food Source Angle")
    print(f"{sun_direction_angle} + {dance_angle_from_vertical} = {food_source_angle}")
    
    food_source_direction_name = [d for d, a in cardinal_angles.items() if a == food_source_angle][0]
    print(f"An angle of {food_source_angle} degrees corresponds to the cardinal direction: {food_source_direction_name}.\n")

    # Step 4: Identify the final location.
    print("Step 4: Identifying the location")
    foraging_location = crop_locations[food_source_direction_name]
    print(f"The field located to the {food_source_direction_name} is the {foraging_location}.")
    print("\nConclusion: The bees are most likely foraging in the strawberry field.")


if __name__ == '__main__':
    solve_bee_foraging_location()