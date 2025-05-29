def average_pooling(matrix, kernel_size):
    if kernel_size == 1:
        # For kernel size 1, return the same matrix formatted to 2 decimal places
        result = [[f"{x:.2f}" for x in row] for row in matrix]
        
        # Print in the required format
        print("<<<")
        for row in result:
            print(" ".join(row))
        print(">>>")

# Input matrix
matrix = [
    [2, 0, 2, 6, 0],
    [5, 1, 9, 1, 3],
    [7, 1, 4, 9, 2],
    [6, 1, 4, 7, 5],
    [1, 4, 0, 5, 3],
    [2, 0, 1, 3, 4]
]

average_pooling(matrix, 1)