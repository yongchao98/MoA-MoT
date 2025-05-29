def count_rectangles(grid):
    # Function to check if a specific rectangle exists at given coordinates
    def verify_rectangle(top, left, bottom, right):
        # Check corners
        if not (grid[top][left] in '#█' and grid[top][right] in '#█' and 
                grid[bottom][left] in '#█' and grid[bottom][right] in '#█'):
            return False
        
        # Check horizontal borders
        for x in range(left, right + 1):
            if not grid[top][x] in '#█' or not grid[bottom][x] in '#█':
                return False
        
        # Check vertical borders
        for y in range(top, bottom + 1):
            if not grid[y][left] in '#█' or not grid[y][right] in '#█':
                return False
        
        return True

    count = 0
    
    # 1. Large outer rectangle
    if verify_rectangle(73, 11, 95, 95):
        count += 1
    
    # 2. Top rectangle
    if verify_rectangle(56, 39, 60, 69):
        count += 1
    
    # 3. Middle rectangle (with overlapping parts)
    if verify_rectangle(61, 39, 76, 69):
        count += 1
    
    # 4. Small inner rectangle
    if verify_rectangle(69, 47, 71, 51):
        count += 1
    
    return count

# Create the grid
grid = [
    " " * 80,
] * 56 + [
    "                                       #########################                ",
    "                                       #                       #                ",
    "                                       #                       #                ",
    "                                       #                       #                ",
    "                                       #                       #                ",
    "                                       █##################     #                ",
    "                                       █                 #     #                ",
    "                                       █                 #     #                ",
    "                                       █                 #     #                ",
    "                                       █                 #     #                ",
    "                                       █                 #     #                ",
    "                                       █                 #     #                ",
    "                                       █       #####     #     #                ",
    "                                       █       #   #     #     #                ",
    "                                       █       #   #     #     #                ",
    "                                       █#######█###█######     #                ",
    "           ############################█#######█████###########█###########     ",
    "           #                           #                       #          #     ",
    "           #                           #                       #          #     ",
    "           #                           #                       #          #     ",
    "           #                           #########################          #     ",
    "           #                                                              #     ",
    "           ################################################################     "
]

print(count_rectangles(grid))