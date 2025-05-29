def find_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    
    # Minimum side length of the square
    min_side = 3
    
    # Iterate over each possible top-left corner of the square
    for i in range(rows - min_side + 1):
        for j in range(cols - min_side + 1):
            # Check for squares of increasing size starting from min_side
            for side in range(min_side, min(rows - i, cols - j) + 1):
                # Check if all characters in the current square are the same
                char = matrix[i][j]
                is_square = True
                for x in range(i, i + side):
                    for y in range(j, j + side):
                        if matrix[x][y] != char:
                            is_square = False
                            break
                    if not is_square:
                        break
                
                # If a valid square is found, record the bottom-right corner
                if is_square:
                    bottom_right = (i + side - 1, j + side - 1)
    
    return bottom_right

# Given character matrix
matrix = [
    ['Q', 'N', 'M', 'K', 'W', 'Y', 'O', 'F'],
    ['B', 'U', 'C', 'Y', 'K', 'F', 'S', 'S'],
    ['I', 'M', 'P', 'M', 'K', 'Q', 'V', 'L'],
    ['Z', 'L', 'U', 'B', 'Q', 'D', 'I', 'Q'],
    ['E', 'Z', 'W', 'H', 'E', 'A', 'V', 'C'],
    ['A', 'P', 'Y', 'Y', 'Y', 'Y', 'A', 'R'],
    ['L', 'F', 'Y', 'Y', 'Y', 'Y', 'U', 'W'],
    ['P', 'F', 'Y', 'Y', 'Y', 'Y', 'G', 'C'],
    ['S', 'M', 'Y', 'Y', 'Y', 'Y', 'H', 'T'],
    ['Q', 'P', 'S', 'Z', 'C', 'K', 'K', 'Y'],
    ['L', 'L', 'O', 'N', 'O', 'G', 'A', 'A'],
    ['W', 'Z', 'G', 'S', 'F', 'Z', 'V', 'A'],
    ['P', 'X', 'P', 'P', 'X', 'N', 'N', 'X'],
    ['Z', 'G', 'Z', 'Q', 'K', 'S', 'Y', 'B']
]

# Find and print the bottom-right corner of the square
bottom_right_corner = find_square(matrix)
print(bottom_right_corner)