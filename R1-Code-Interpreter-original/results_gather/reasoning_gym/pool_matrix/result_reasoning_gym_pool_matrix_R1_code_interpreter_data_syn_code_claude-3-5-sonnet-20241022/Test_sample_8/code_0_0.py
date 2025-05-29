def average_pooling(matrix, kernel_size):
    if not matrix or not matrix[0]:
        return []
    
    rows = len(matrix)
    cols = len(matrix[0])
    
    # For kernel size 1, output will be same as input
    if kernel_size == 1:
        # Format each number to 2 decimal places
        result = [[round(x, 2) for x in row] for row in matrix]
        
        # Print in the required format
        print("<<<")
        for row in result:
            print(" ".join(f"{x:.2f}" for x in row))
        print(">>>")

# Input matrix
matrix = [
    [6, 7, 5, 8, 8, 2, 1, 7, 2, 3],
    [7, 3, 9, 0, 1, 5, 9, 8, 4, 4],
    [8, 6, 3, 3, 7, 7, 1, 3, 5, 2],
    [2, 3, 7, 9, 3, 2, 2, 5, 0, 1],
    [7, 9, 3, 8, 5, 9, 3, 2, 6, 9],
    [9, 3, 5, 4, 5, 5, 8, 0, 1, 1],
    [4, 8, 9, 2, 5, 7, 8, 6, 4, 2],
    [9, 7, 3, 9, 0, 5, 5, 5, 0, 3],
    [2, 4, 9, 7, 6, 7, 7, 3, 4, 9]
]

average_pooling(matrix, 1)