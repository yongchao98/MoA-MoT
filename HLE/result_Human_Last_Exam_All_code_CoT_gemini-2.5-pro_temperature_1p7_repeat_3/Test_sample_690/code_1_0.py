import math

def flawed_leftmost_search(row_slice):
    """
    Implements the flawed binary search for the leftmost non-zero element.
    Returns the index in the slice, or -1 if not found.
    """
    left, right = 0, len(row_slice) - 1
    leftmost_col = -1
    while left <= right:
        mid = (left + right) // 2
        if row_slice[mid] == 0:
            left = mid + 1
        else:
            leftmost_col = mid
            right = mid - 1
    return leftmost_col

def flawed_rightmost_search(row_slice):
    """
    Implements the flawed binary search for the rightmost non-zero element.
    Returns the index in the slice, or -1 if not found.
    """
    left, right = 0, len(row_slice) - 1
    rightmost_col = -1
    while left <= right:
        mid = (left + right) // 2
        if row_slice[mid] == 0:
            right = mid - 1
        else:
            rightmost_col = mid
            left = mid + 1
    return rightmost_col

def calculate_bandwidth_flawed(matrix):
    """
    Implements the flawed algorithm described in the problem.
    Matrix indices are 0-based.
    """
    n = len(matrix)
    if n == 0:
        return 0
    
    overall_bandwidth = 0
    for i in range(n):
        # Step 3.b: Find leftmost non-zero in A[i, 0...i]
        left_slice = matrix[i][:i+1]
        # flawed_leftmost_search returns index within the slice
        leftmost_idx_in_slice = flawed_leftmost_search(left_slice)
        # if no non-zero is found, default to the diagonal
        leftmost_j = leftmost_idx_in_slice if leftmost_idx_in_slice != -1 else i

        # Step 3.c: Find rightmost non-zero in A[i, i...n-1]
        right_slice = matrix[i][i:]
        rightmost_idx_in_slice = flawed_rightmost_search(right_slice)
        # if no non-zero is found, default to the diagonal
        # Index needs to be offset by `i` to get the original matrix column index
        rightmost_j = (rightmost_idx_in_slice + i) if rightmost_idx_in_slice != -1 else i

        # Step 3.d: Calculate distance and row bandwidth
        dist_left = i - leftmost_j
        dist_right = rightmost_j - i
        row_bandwidth = max(dist_left, dist_right)

        # Step 3.e: Update overall bandwidth
        if row_bandwidth > overall_bandwidth:
            overall_bandwidth = row_bandwidth
            
    return overall_bandwidth

def calculate_bandwidth_correct(matrix):
    """A correct O(n^2) brute-force algorithm for comparison."""
    n = len(matrix)
    max_dist = 0
    true_bw_i, true_bw_j = -1, -1
    for i in range(n):
        for j in range(n):
            if matrix[i][j] != 0:
                dist = abs(i - j)
                if dist > max_dist:
                    max_dist = dist
                    true_bw_i, true_bw_j = i, j
    return max_dist, true_bw_i, true_bw_j

# A simple non-symmetric matrix that demonstrates the failure.
# The binary search will fail to find the non-zero element A[0][2].
test_matrix = [
    [0, 0, 5, 0],
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1]
]

# Run both algorithms
flawed_result = calculate_bandwidth_flawed(test_matrix)
correct_result, i, j = calculate_bandwidth_correct(test_matrix)

# Print the results and explanation
print("Demonstrating the algorithm's flaw:")
print(f"Test Matrix:\n{test_matrix[0]}\n{test_matrix[1]}\n{test_matrix[2]}\n{test_matrix[3]}\n")

print("The flawed algorithm's binary search fails to find the rightmost non-zero element (5 at index 2) in row 0, `[0, 0, 5, 0]`.")
print(f"The flawed algorithm calculates a bandwidth of: {flawed_result}")

print("\nThe correct (brute-force) algorithm finds:")
print(f"The true bandwidth is determined by the element A[{i}][{j}].")
print(f"|{i} - {j}| = {correct_result}")
print(f"The true bandwidth is: {correct_result}")
