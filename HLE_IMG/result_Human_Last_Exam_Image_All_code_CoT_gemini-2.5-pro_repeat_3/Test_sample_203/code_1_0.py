def solve_heesch_number():
    """
    This function determines and prints the Heesch number for the three polygons.

    Based on geometric analysis:
    1. The three polygons are identical.
    2. The polygon shape has four outward protrusions and no corresponding indentations.
    3. This means an identical copy cannot be placed snugly against any of the protruded sides without overlapping.
    4. Therefore, a complete "first layer" of surrounding polygons cannot be formed.
    5. By definition, if the first layer cannot be completed, the Heesch number is 0.
    6. Since all three polygons are the same, the Heesch number for each is 0.
    """
    
    # Heesch number for the first polygon
    heesch_number_1 = 0
    
    # Heesch number for the second polygon
    heesch_number_2 = 0
    
    # Heesch number for the third polygon
    heesch_number_3 = 0
    
    # The problem asks for the answers separated by commas.
    print(f"{heesch_number_1},{heesch_number_2},{heesch_number_3}")

solve_heesch_number()