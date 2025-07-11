def solve_snail_puzzle():
    """
    Solves the snail puzzle by calculating the maximal distance.

    The problem can be solved by considering three parts of the snail's journey,
    which are structured to maximize the distance observed.
    
    1. First Part: A 3-meter journey observed over the first part of the timeline.
    2. Second Part: A 4-meter journey observed cleverly in the middle.
    3. Third Part: A 3-meter journey observed over the final part of the timeline.

    The total maximal distance is the sum of these parts.
    """
    
    # Distance from the first part of the journey.
    part1_distance = 3
    
    # Distance from the clever middle part of the journey.
    part2_distance = 4
    
    # Distance from the final part of the journey.
    part3_distance = 3
    
    # The maximal total distance is the sum of these parts.
    total_distance = part1_distance + part2_distance + part3_distance
    
    # Output the final calculation clearly.
    print("The maximal distance is the sum of three conceptual parts of the journey.")
    print(f"Part 1: {part1_distance} meters")
    print(f"Part 2: {part2_distance} meters")
    print(f"Part 3: {part3_distance} meters")
    print("The final equation is:")
    print(f"{part1_distance} + {part2_distance} + {part3_distance} = {total_distance}")

solve_snail_puzzle()