def rotate_540(matrix):
    # For 540 degrees (180 degrees), we just need to reverse the matrix
    # both vertically and horizontally
    n = len(matrix)
    
    # Create a string representation of the rotated matrix
    result = []
    for i in range(n-1, -1, -1):  # reverse rows
        row = matrix[i][::-1]  # reverse elements in each row
        result.append(' '.join(map(str, row)))
    
    return '<<<' + '\n'.join(result) + '>>>'

# Input matrix
matrix = [
    [7, 5, 9, 2],
    [3, 3, 7, 2],
    [1, 5, 0, 3],
    [1, 7, 4, 2]
]

print(rotate_540(matrix))