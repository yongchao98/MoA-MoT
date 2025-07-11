import math

def calculate_path_danger(path_coords, dragon_coord):
    """
    Calculates the danger of a path.
    The danger is inversely related to the minimum distance to the dragon.
    A smaller minimum distance means a more dangerous path.
    """
    min_dist = float('inf')
    closest_point = None

    # Using Manhattan distance for grid-based movement
    for point in path_coords:
        distance = abs(point[0] - dragon_coord[0]) + abs(point[1] - dragon_coord[1])
        if distance < min_dist:
            min_dist = distance
            closest_point = point
            
    return min_dist, closest_point

def main():
    """
    Main function to analyze the map and find the least dangerous path.
    """
    # Coordinates are in (row, column) format
    start_pos = (3, 25)
    gold_pos = (3, 4)
    dragon_pos = (10, 25)

    print("Analyzing possible paths...")
    print(f"Adventurer Start: {start_pos}, Gold: {gold_pos}, Dragon: {dragon_pos}\n")

    # --- Path 1: The High Road (Straight Left) ---
    # This path moves left from the start to the gold along row 3.
    # It passes through an unlit hallway (#), which is considered traversable.
    path1_coords = [(3, c) for c in range(start_pos[1], gold_pos[1] - 1, -1)]
    path1_moves = "L"

    # --- Path 2: The Low Road (Down, Up, Left) ---
    # This path goes down through the door at (4, 25), into the vertical corridor.
    # We assume a path through the corridor like: (5,24) -> Up to (3,24) -> Left to gold.
    # This is a more complex path that brings the adventurer closer to the dragon.
    # Initial segment: Down from (3,25) to (4,25), then into the corridor at (4,24).
    path2_initial_segment = [(3, 25), (4, 25), (4, 24)] 
    # The rest of the path would go up and left.
    path2_moves = "DUL"

    # Note: A detailed analysis shows the Low Road (Path 2) is a likely dead-end
    # because of the wall at row 4. However, for the sake of a comprehensive danger
    # analysis, we will calculate its danger as if it were possible.

    # --- Danger Calculation ---
    min_dist1, closest_pt1 = calculate_path_danger(path1_coords, dragon_pos)
    # For path 2, we only need to check its most dangerous part, the initial segment.
    min_dist2, closest_pt2 = calculate_path_danger(path2_initial_segment, dragon_pos)

    print("--- Path 1 (High Road) Analysis ---")
    print(f"Move Sequence: {path1_moves}")
    print(f"Closest point to dragon: {closest_pt1}")
    print(f"The minimum distance to the dragon is calculated as: |{closest_pt1[0]} - {dragon_pos[0]}| + |{closest_pt1[1]} - {dragon_pos[1]}| = {min_dist1}")
    print(f"Minimum Distance: {min_dist1}\n")


    print("--- Path 2 (Low Road) Analysis ---")
    print(f"Move Sequence: {path2_moves}")
    print(f"Closest point to dragon: {closest_pt2}")
    print(f"The minimum distance to the dragon is calculated as: |{closest_pt2[0]} - {dragon_pos[0]}| + |{closest_pt2[1]} - {dragon_pos[1]}| = {min_dist2}")
    print(f"Minimum Distance: {min_dist2}\n")

    # --- Conclusion ---
    print("--- Conclusion ---")
    if min_dist1 > min_dist2:
        print(f"Path 1 is less dangerous because its minimum distance to the dragon ({min_dist1}) is greater than Path 2's ({min_dist2}).")
        final_path = path1_moves
    else:
        print(f"Path 2 is less dangerous because its minimum distance to the dragon ({min_dist2}) is greater than Path 1's ({min_dist1}).")
        final_path = path2_moves
        
    print(f"\nThe least dangerous path is: {final_path}")
    print(f"<<<{final_path}>>>")

if __name__ == "__main__":
    main()