def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = None

    # Check for squares starting from each position
    for i in range(rows):
        for j in range(cols):
            # Try to form a square with top-left corner at (i, j)
            current_char = matrix[i][j]
            size = 2  # Start checking from size 3x3
            while i + size < rows and j + size < cols:
                # Check if all characters in the square are the same
                valid_square = True
                for x in range(i, i + size + 1):
                    for y in range(j, j + size + 1):
                        if matrix[x][y] != current_char:
                            valid_square = False
                            break
                    if not valid_square:
                        break
                
                if valid_square:
                    if size + 1 > max_size:
                        max_size = size + 1
                        bottom_right_corner = (i + size, j + size)
                else:
                    break
                size += 1

    return bottom_right_corner

matrix = [
    ['L', 'H', 'P', 'T', 'T', 'K', 'A', 'S'],
    ['B', 'W', 'W', 'W', 'D', 'I', 'J', 'W'],
    ['U', 'W', 'W', 'W', 'A', 'A', 'E', 'A'],
    ['O', 'W', 'W', 'W', 'T', 'T', 'Y', 'S'],
    ['D', 'M', 'Z', 'V', 'D', 'Z', 'Y', 'V'],
    ['S', 'C', 'T', 'W', 'H', 'C', 'Z', 'P'],
    ['K', 'M', 'N', 'Q', 'E', 'H', 'D', 'N'],
    ['E', 'K', 'S', 'X', 'A', 'U', 'I', 'D'],
    ['F', 'F', 'D', 'S', 'V', 'J', 'P', 'J'],
    ['I', 'I', 'S', 'I', 'T', 'I', 'V', 'Z'],
    ['R', 'Z', 'L', 'W', 'Y', 'J', 'D', 'E']
]

result = find_largest_square(matrix)
print(result)