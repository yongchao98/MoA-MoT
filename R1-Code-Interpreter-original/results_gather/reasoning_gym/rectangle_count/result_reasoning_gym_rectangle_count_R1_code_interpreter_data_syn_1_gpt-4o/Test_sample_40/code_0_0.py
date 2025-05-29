# Define the grid as a list of strings
grid = [
    "                                       #########################################",
    "                                       #                                       #",
    "                                       #                                       #",
    "                                       #                                       #",
    "                                       #                                       #",
    "                                       #                                       #",
    "                                       #                                       #",
    "                                       #                                       #",
    "                                       #                                       #",
    "                                       #                                       #",
    "                                       #                                       #",
    "                                       #                                       #",
    "                                       #                                       #",
    "                                       #                                       #",
    "                                       #                                       #",
    "                                       #                                       #",
    "                                       #                                       #",
    "                                       #                                       #",
    "                                       #                                       #",
    "                                       #                                       #",
    "                                       #                                       #",
    "                                       #                                       #",
    "                                       #                                       #",
    "                                       #                                       #",
    "                                       #                                       #",
    "                                       #                                       #",
    "                                       #                                       #",
    "                                       #                                       #",
    "                                       #                                    ####",
    "                                       #                                    # ##",
    "                                       #                                    # ##",
    "                                       #                                    # ##",
    "                                       #                                    # ##",
    "                                       #                                    # ##",
    "                                       #                                    # ##",
    "                                       #                                    # ##",
    "                                       #                                    # ##",
    "                                       #                                    # ##",
    "                                       #                                    # ##",
    "                                       #                                    # ##",
    "                                       #                                    # ##",
    "                                       #                                    # ##",
    "                                       #                                    # ##",
    "                                       #                                    # ##",
    "                                       #                                    # ##",
    "                                       #####################################█#█#",
    "                                                                            # # ",
    "                                                                            # # ",
    "                                                                            # # ",
    "                                                                            # # ",
    "                                                                            # # ",
    "                                                                            # # ",
    "                                                                            ### "
]

# Initialize variables
rectangles = 0
overlap_detected = False

# Scan the grid
for row in grid:
    if '█' in row:
        overlap_detected = True
        break

# Count rectangles
if overlap_detected:
    rectangles = 2
else:
    rectangles = 1

# Output the result
print(rectangles)