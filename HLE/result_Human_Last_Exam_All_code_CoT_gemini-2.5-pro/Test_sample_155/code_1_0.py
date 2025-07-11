def solve_and_decode():
    """
    Parses the visual encoding and prints the decoded word shape.
    The word is found in the "negative space" of the original image.
    """
    verbatim_text = """t b b t t t t b b b
b
b
b
t
b
t b t t b b b b b b
                  b
                  b
                  t
                  t
                  b
b b t t t t b t b b"""

    lines = verbatim_text.split('\n')
    
    # Determine the maximum width to create a consistent grid
    max_width = 0
    for line in lines:
        if len(line) > max_width:
            max_width = len(line)

    # Create a grid, padding all lines to the same width
    grid = []
    for line in lines:
        grid.append(list(line.ljust(max_width)))

    # Create a condensed matrix using only the even-indexed columns,
    # which contain the actual characters of the drawing.
    condensed_matrix = []
    height = len(grid)
    # The width of the condensed matrix is half the original width.
    condensed_width = (max_width + 1) // 2 
    
    for r in range(height):
        matrix_row = []
        for c in range(condensed_width):
            # Map matrix column 'c' to original grid column 'c*2'
            original_col_index = c * 2
            if original_col_index < len(grid[r]):
                matrix_row.append(grid[r][original_col_index])
            else:
                matrix_row.append(' ') # Should not happen with ljust
        condensed_matrix.append(matrix_row)
        
    print("The decoded message is formed by the shapes below:")
    # Print the inverted matrix to reveal the word in the negative space.
    for r in range(height):
        output_line = ""
        for c in range(condensed_width):
            char = condensed_matrix[r][c]
            if char in ('t', 'b'):
                # Replace original characters with space
                output_line += " "
            else:
                # Replace original spaces with '#' to show the word
                output_line += "#"
        print(output_line)

solve_and_decode()