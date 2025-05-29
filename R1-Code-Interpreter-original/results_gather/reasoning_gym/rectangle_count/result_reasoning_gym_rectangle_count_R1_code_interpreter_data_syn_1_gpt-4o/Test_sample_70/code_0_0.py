# Define the grid as a list of strings
grid = [
    "                                                 ############                   ",
    "                                                 #          #                   ",
    "                                                 #          #                   ",
    "                                                 #          #                   ",
    "                                                 #          #                   ",
    "                                                 #          #                   ",
    "                                                 #          #                   ",
    "                                                 #          #                   ",
    "                                                 #          #                   ",
    "                                                 #          #                   ",
    "                                                 #          #                   ",
    "                                                 #          #                   ",
    "                                                 #          #                   ",
    "                                                 #          #             ###   ",
    "                                                 #          #             # #   ",
    "                                                 #          #             # #   ",
    "                                             ####█########  #             # #   ",
    "                                             #   #       #  #             # #   ",
    "                                             #   #       #  #             # #   ",
    "                                             #   #       #  #             # #   ",
    "                                             ####█████████###             # #   ",
    "                                                                          # #   ",
    "                                                                          # #   ",
    "                                                                          ###   "
]

def count_rectangles(grid):
    # Initialize variables
    rectangles = 0
    in_rectangle = False
    overlap_detected = False

    # Iterate over each line in the grid
    for line in grid:
        # Check for the start of a rectangle
        if '#' in line or '█' in line:
            if not in_rectangle:
                rectangles += 1
                in_rectangle = True
            # Check for overlap
            if '█' in line:
                overlap_detected = True
        else:
            in_rectangle = False

    # If overlap was detected, add one more rectangle
    if overlap_detected:
        rectangles += 1

    return rectangles

# Count the rectangles
total_rectangles = count_rectangles(grid)
print(total_rectangles)