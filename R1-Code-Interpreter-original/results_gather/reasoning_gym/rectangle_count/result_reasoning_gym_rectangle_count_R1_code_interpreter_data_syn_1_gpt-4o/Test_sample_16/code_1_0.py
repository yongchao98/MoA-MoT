def count_rectangles(grid):
    # Function to parse the grid and count rectangles
    def is_rectangle_char(c):
        return c == '#' or c == '█'

    def find_rectangles_in_line(line):
        rectangles = []
        in_rectangle = False
        start = 0
        for i, char in enumerate(line):
            if is_rectangle_char(char) and not in_rectangle:
                in_rectangle = True
                start = i
            elif not is_rectangle_char(char) and in_rectangle:
                in_rectangle = False
                rectangles.append((start, i - 1))
        if in_rectangle:
            rectangles.append((start, len(line) - 1))
        return rectangles

    # Parse the grid
    rectangles = []
    for y, line in enumerate(grid):
        line_rectangles = find_rectangles_in_line(line)
        for start, end in line_rectangles:
            rectangles.append((start, y, end, y))

    # Merge rectangles vertically
    merged_rectangles = []
    while rectangles:
        rect = rectangles.pop(0)
        start_x, start_y, end_x, end_y = rect
        merged = False
        for i, other in enumerate(rectangles):
            if other[0] == start_x and other[2] == end_x and other[1] == end_y + 1:
                end_y = other[3]
                rectangles.pop(i)
                rectangles.append((start_x, start_y, end_x, end_y))
                merged = True
                break
        if not merged:
            merged_rectangles.append((start_x, start_y, end_x, end_y))

    return len(merged_rectangles)

# Define the grid
grid = [
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                            #########           ",
    "                                                            #       #           ",
    "                                                            #       #           ",
    "                                                            #       #           ",
    "                    ###########                             #       #           ",
    "                    #         #                             #       #           ",
    "                    #         #                             #       #           ",
    "                    #         #                             #       #           ",
    "                    #         #                             #       #           ",
    "                    #         #                             #       #           ",
    "                    #         #                             #       #           ",
    "                    #         #                             #########           ",
    "                    #         #                                                 ",
    "                    #         #                       #####                     ",
    "                    #         #                       #   #                     ",
    "                    #         #                       #   #                     ",
    "                    #         #                       ####█####                 ",
    "                    #         #                       ##  #   #                 ",
    "                    #         #                       ##  #   #                 ",
    "                    #         #                       ##  #   #                 ",
    "                    #         #                       ##  #   #                 ",
    "                    #         #                       ##  #   #                 ",
    "                    #         #                       ##  #   #                 ",
    "                    #         #                       ##  #   #                 ",
    "                    #         #                       ##  #   #                 ",
    "             #######█#########█#######################██##█###█################ ",
    "             #      #         #                       ##  #   #               # ",
    "             #      #         #                       ##  #   #               # ",
    "             #      #         #                       ##  #   #               # ",
    "             #      #         #                       ##  #   #               # ",
    "             #      #         #                       ##  #   #               # ",
    "             #      #         #                       ##  #   #               # ",
    "             #      #         #                       ##  #   #               # ",
    " ############█######█#########█#######################██##█###█############## # ",
    " #           #      #         #                       ##  #   #             # # ",
    " #           #      #         #                       ##  #   #             # # ",
    " #           #      #         #                       ##  #   #             # # ",
    " #           #      #         #                       ##  #   #             # # ",
    " #           #      #         #                       ##  #   #             # # ",
    " #           #      #         #                       ##  #   #             # # ",
    " #           #######█#########█#######################██##█###█#############█## ",
    " #                  #         #                       ##  #   #             #   ",
    " #                  #         #                       ##  #   #             #   ",
    " #                  #         #                       ##  #   #             #   ",
    " #                  #         #                       ##  #   #             #   ",
    " #                  #         #                       ##  #   #             #   ",
    " #                  #         #                       #█###   #             #   ",
    " #                  #         #                        #      #             #   ",
    " #                  #         #                        #      #             #   ",
    " #                  #         #                        #      #             #   ",
    " #                  #         #                        #      #             #   ",
    " #                  #         #                        #      #             #   ",
    " #                  #         #                        #      #             #   ",
    " #                  #         #                        #      #             #   ",
    " #                  #         #                        #      #             #   ",
    " #                  #         #                        #      #             #   ",
    " #                  #         #                        #      #             #   ",
    " #                  #         #                        #      #             #   ",
    " #                  #         #                        #      #             #   ",
    " #                  #         #                        #      #             #   ",
    " #                  #         #                        #      #             #   ",
    " #                  #         #                        #      #             #   ",
    " #                  #         #                        #      #     ####### #   ",
    " #                  #         #                        #      #  ###█#####█##   ",
    " #                  #         #                        #      #  #  #     ###   ",
    " #                  #         #                        #      #  #  #     ###   ",
    " #                  #         #                        #      #  #  #     ###   ",
    " #                  #         #                        #      #  ###███████##   ",
    " #                  #         #                        #      #             #   ",
    " #                  #         #         ###########    #      #             #   ",
    " #                  #         #         #         #    ########             #   ",
    " #                  #         #         #         #                         #   ",
    " #                  #         #         #         #                         #   ",
    " #                  #         #         #         #                         #   ",
    " #                  #         #         #         #                         #   ",
    " #                  #         #         #         #                         #   ",
    " ###################█#########█#########███████████##########################   ",
    "                    ###########                                                 "
]

# Count rectangles
rectangle_count = count_rectangles(grid)
print(rectangle_count)