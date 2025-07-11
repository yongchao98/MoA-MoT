import math

def flawed_bandwidth_algorithm(matrix):
    """
    Implements the flawed algorithm described in the prompt.
    Uses 0-based indexing.
    """
    n = len(matrix)
    overall_bandwidth = 0

    print("Running the flawed algorithm from the prompt...")

    for i in range(n):
        # Step 3b: Flawed binary search for the leftmost non-zero in A[i, 0...i]
        leftmost_col = -1
        low, high = 0, i
        current_best = i + 1
        
        while low <= high:
            mid = (low + high) // 2
            if matrix[i][mid] != 0:
                current_best = mid
                high = mid - 1
            else:
                low = mid + 1
        leftmost_col = current_best if current_best <= i else -1

        # Step 3c: Flawed binary search for the rightmost non-zero in A[i, i...n-1]
        rightmost_col = -1
        low, high = i, n - 1
        current_best = i - 1
        
        while low <= high:
            mid = (low + high) // 2
            if matrix[i][mid] != 0:
                current_best = mid
                low = mid + 1
            else:
                high = mid - 1
        rightmost_col = current_best if current_best >= i else -1

        # Step 3d: Calculate row bandwidth
        dist_left = 0
        if leftmost_col != -1 and leftmost_col <= i:
            dist_left = i - leftmost_col

        dist_right = 0
        if rightmost_col != -1 and rightmost_col >= i:
            dist_right = rightmost_col - i
        
        row_bandwidth = max(dist_left, dist_right)
        
        # Step 3e: Update overall bandwidth
        if row_bandwidth > overall_bandwidth:
            overall_bandwidth = row_bandwidth
            
        print(f"Row {i}: Leftmost found at {leftmost_col}, Rightmost found at {rightmost_col}. "
              f"Distances: ({i}-{leftmost_col}={dist_left}, {rightmost_col}-{i}={dist_right}). "
              f"Row bandwidth = {row_bandwidth}. Overall bandwidth = {overall_bandwidth}")

    return overall_bandwidth

def correct_bandwidth_algorithm(matrix):
    """
    Correctly computes the bandwidth using a simple O(n^2) brute-force search.
    """
    n = len(matrix)
    bandwidth = 0
    for i in range(n):
        for j in range(n):
            if matrix[i][j] != 0:
                distance = abs(i - j)
                if distance > bandwidth:
                    bandwidth = distance
    return bandwidth

# A symmetric matrix where the flawed algorithm will fail.
# The non-zero elements at (0,4) and (4,0) are separated by zeros.
test_matrix = [
    [1, 0, 0, 0, 5],
    [0, 1, 0, 0, 0],
    [0, 0, 1, 0, 0],
    [0, 0, 0, 1, 0],
    [5, 0, 0, 0, 1]
]

print("Matrix A:")
for row in test_matrix:
    print(row)
print("-" * 20)

# Run the flawed algorithm
flawed_result = flawed_bandwidth_algorithm(test_matrix)
print("-" * 20)
print(f"Final bandwidth from flawed algorithm: {flawed_result}\n")

# Run the correct algorithm
correct_result = correct_bandwidth_algorithm(test_matrix)
print(f"Correct bandwidth (from brute-force): {correct_result}")
print("-" * 20)

# Final conclusion based on the analysis
print("Analysis Summary:")
print("1. Correctness: The algorithm is INCORRECT because binary search fails on non-monotonic data (rows with interleaved zeros and non-zeros).")
print("2. Complexity: The described procedure involves a loop of size n with O(log n) work inside, giving a total complexity of O(n*log(n)).")
print("This matches option C.")

<<<C>>>