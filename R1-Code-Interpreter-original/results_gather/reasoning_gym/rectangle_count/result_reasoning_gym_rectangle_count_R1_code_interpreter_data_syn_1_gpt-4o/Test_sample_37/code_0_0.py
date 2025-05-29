# Define the grid as a list of strings
grid = [
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "            ########################                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            ########################                                            ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
]

# Initialize variables
rectangle_count = 0
inside_rectangle = False

# Iterate over each line in the grid
for line in grid:
    # Check if the line contains a rectangle border
    if '#' in line:
        if not inside_rectangle:
            # We found the top border of a new rectangle
            rectangle_count += 1
            inside_rectangle = True
    else:
        # If we encounter a line without '#', we are outside a rectangle
        inside_rectangle = False

# Output the total count of rectangles
print(rectangle_count)