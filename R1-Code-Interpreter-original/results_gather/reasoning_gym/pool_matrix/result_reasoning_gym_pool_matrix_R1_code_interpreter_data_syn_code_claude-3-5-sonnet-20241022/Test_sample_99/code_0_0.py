import numpy as np

# Input matrix
matrix = [
    [3, 7, 9, 9, 5, 0, 2],
    [8, 1, 8, 0, 5, 4, 0],
    [0, 1, 1, 2, 7, 5, 8],
    [7, 5, 0, 7, 0, 3, 8],
    [8, 8, 6, 8, 3, 2, 4],
    [5, 3, 5, 3, 8, 0, 9],
    [1, 8, 2, 2, 6, 5, 6]
]

kernel_size = 3

def average_pooling(matrix, kernel_size):
    n = len(matrix)
    output_size = n // kernel_size
    result = []
    
    for i in range(output_size):
        row = []
        for j in range(output_size):
            # Extract the kernel region
            values = []
            for ki in range(kernel_size):
                for kj in range(kernel_size):
                    values.append(matrix[i*kernel_size + ki][j*kernel_size + kj])
            # Calculate average and round to 2 decimal places
            avg = round(sum(values) / (kernel_size * kernel_size), 2)
            row.append(avg)
        result.append(row)
    
    # Format the output
    output = ""
    for row in result:
        output += " ".join(str(x) for x in row) + "\n"
    print(output.strip())

average_pooling(matrix, kernel_size)