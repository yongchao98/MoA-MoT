import numpy as np

def max_pooling(matrix, window_size):
    n = len(matrix)
    result_size = n - window_size + 1
    result = np.zeros((result_size, result_size), dtype=int)
    
    for i in range(result_size):
        for j in range(result_size):
            window = matrix[i:i+window_size, j:j+window_size]
            result[i, j] = np.max(window)
    
    return result

matrix = np.array([
    [6, 3, 8, 4, 20],
    [-8, -8, 0, 19, -10],
    [-7, 15, 2, -2, 10],
    [9, -9, 10, 7, 13],
    [-5, 11, -9, -1, 12]
])

result = max_pooling(matrix, 3)
print(result)