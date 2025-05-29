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
    [-2, 12, -3],
    [11, 20, -3],
    [-9, -4, 10]
])

window_size = 2
result = max_pooling(matrix, window_size)
print(result)