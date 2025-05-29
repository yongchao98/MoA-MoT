# Input matrix
matrix = [
    [3, 2, 3, 6, 2, 0, 1, 1, 3],
    [3, 9, 5, 4, 0, 1, 9, 7, 1],
    [3, 4, 8, 8, 9, 3, 4, 0, 0],
    [7, 3, 1, 3, 0, 4, 6, 6, 2]
]

def max_pooling(matrix, kernel_size):
    rows = len(matrix)
    cols = len(matrix[0])
    
    # Calculate output dimensions
    output_rows = rows // kernel_size
    output_cols = cols // kernel_size
    
    result = []
    
    # Process each region
    for i in range(output_rows):
        row = []
        for j in range(output_cols):
            # Get the maximum value in the current kernel region
            current_max = float('-inf')
            for ki in range(kernel_size):
                for kj in range(kernel_size):
                    r = i * kernel_size + ki
                    c = j * kernel_size + kj
                    if r < rows and c < cols:
                        current_max = max(current_max, matrix[r][c])
            row.append(current_max)
        result.append(row)
    
    # Format output
    formatted_result = [[f"{x:.2f}" for x in row] for row in result]
    
    # Print in required format
    print("<<<", end="")
    for i, row in enumerate(formatted_result):
        print(" ".join(row), end="")
        if i < len(formatted_result) - 1:
            print()
    print(">>>")

# Perform max pooling with kernel size 3
max_pooling(matrix, 3)