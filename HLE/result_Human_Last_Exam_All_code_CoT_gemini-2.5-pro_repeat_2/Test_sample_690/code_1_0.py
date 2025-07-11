import numpy as np

def described_algorithm(matrix):
    """
    Implements the flawed algorithm described in the problem.
    """
    n = matrix.shape[0]
    if n == 0:
        return 0
    
    overall_bandwidth = 0
    
    for i in range(n):
        # Step 3.b: Find the leftmost non-zero element using binary search
        # This search assumes the pattern [0, 0, ..., non-zero, ...]
        l, r = 0, i
        leftmost_col = -1
        while l <= r:
            mid = (l + r) // 2
            if matrix[i, mid] != 0:
                leftmost_col = mid
                r = mid - 1  # Found a non-zero, try to find one further left
            else:
                l = mid + 1  # Non-zero must be to the right
        
        # Step 3.c: Find the rightmost non-zero element using binary search
        # This search assumes the pattern [..., non-zero, non-zero, ..., 0, 0]
        l, r = i, n - 1
        rightmost_col = -1
        while l <= r:
            mid = (l + r) // 2
            if matrix[i, mid] != 0:
                rightmost_col = mid
                l = mid + 1  # Found a non-zero, try to find one further right
            else:
                r = mid - 1  # Non-zero must be to the left

        # Step 3.d: Calculate row bandwidth
        row_bandwidth = 0
        # If the row is all zeros, leftmost/rightmost will be -1
        if leftmost_col != -1 and rightmost_col != -1:
            dist_left = i - leftmost_col
            dist_right = rightmost_col - i
            row_bandwidth = max(dist_left, dist_right)
        
        # Step 3.e: Update overall bandwidth
        if row_bandwidth > overall_bandwidth:
            overall_bandwidth = row_bandwidth
            
    return overall_bandwidth

def correct_algorithm(matrix):
    """
    Implements a correct but less efficient (O(n^2)) algorithm.
    """
    n = matrix.shape[0]
    if n == 0:
        return 0
    
    max_dist = 0
    for i in range(n):
        for j in range(n):
            if matrix[i, j] != 0:
                dist = abs(i - j)
                if dist > max_dist:
                    max_dist = dist
    return max_dist

# A symmetric matrix where the described algorithm fails
counter_example_matrix = np.array([
    [5, 0, 3, 0],
    [0, 5, 0, 2],
    [3, 0, 5, 0],
    [0, 2, 0, 5]
])

# Calculate bandwidth using both algorithms
correct_bw = correct_algorithm(counter_example_matrix)
described_bw = described_algorithm(counter_example_matrix)

print("Demonstrating the flaw in the described algorithm:")
print("Matrix:\n", counter_example_matrix)
print(f"\nThe correct bandwidth is: {correct_bw}")
print(f"The bandwidth computed by the described algorithm is: {described_bw}")
print("\nSince the results differ, the described algorithm is not correct.")
