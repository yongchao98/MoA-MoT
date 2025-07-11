import math

def find_leftmost_nonzero(row_data, start, end):
    """
    Binary search for the leftmost non-zero element in a slice of a row.
    Searches in row_data[start...end].
    """
    left, right = start, end
    leftmost_col = -1

    while left <= right:
        mid = (left + right) // 2
        if row_data[mid] != 0:
            leftmost_col = mid
            right = mid - 1  # Try to find a non-zero element further to the left
        else:
            left = mid + 1   # Non-zero element must be to the right
            
    return leftmost_col

def find_rightmost_nonzero(row_data, start, end):
    """
    Binary search for the rightmost non-zero element in a slice of a row.
    Searches in row_data[start...end].
    """
    left, right = start, end
    rightmost_col = -1

    while left <= right:
        mid = (left + right) // 2
        if row_data[mid] != 0:
            rightmost_col = mid
            left = mid + 1   # Try to find a non-zero element further to the right
        else:
            right = mid - 1  # Non-zero element must be to the left

    return rightmost_col

def calculate_bandwidth_and_show_steps(matrix):
    """
    Calculates the bandwidth of a square matrix using the described algorithm
    and prints the intermediate steps for clarity.
    """
    n = len(matrix)
    if n == 0:
        return 0
    
    overall_bandwidth = 0
    print("Analyzing matrix:")
    for row in matrix:
        print(f"  {row}")
    print("-" * 30)

    # 2. Iterate through each row of the matrix
    for i in range(n):
        row = matrix[i]
        
        # 3b. Find the leftmost non-zero element
        leftmost_col = find_leftmost_nonzero(row, 0, i)
        
        # 3c. Find the rightmost non-zero element
        rightmost_col = find_rightmost_nonzero(row, i, n - 1)
        
        # 3d. Calculate the distance from the diagonal
        dist_left = 0
        if leftmost_col != -1:
            dist_left = i - leftmost_col
            
        dist_right = 0
        if rightmost_col != -1:
            dist_right = rightmost_col - i
            
        row_bandwidth = max(dist_left, dist_right)
        
        print(f"Row {i}:")
        print(f"  Left search in columns [0, {i}]: found non-zero at column {leftmost_col if leftmost_col != -1 else 'None'}")
        print(f"  Right search in columns [{i}, {n-1}]: found non-zero at column {rightmost_col if rightmost_col != -1 else 'None'}")
        print(f"  Row distances from diagonal: left = {dist_left}, right = {dist_right}")
        print(f"  Row bandwidth = max({dist_left}, {dist_right}) = {row_bandwidth}")

        # 3e. Update the overall bandwidth
        if row_bandwidth > overall_bandwidth:
            overall_bandwidth = row_bandwidth
            print(f"  New overall bandwidth updated to {overall_bandwidth}")
        print("-" * 30)
        
    print(f"Final Result:")
    print(f"The bandwidth of the matrix is the maximum row bandwidth found, which is {overall_bandwidth}.")
    
# --- User Execution Area ---
# A sample symmetric matrix to test the algorithm
# The true bandwidth is determined by A[1, 4], where |1 - 4| = 3.
sample_matrix = [
    [5, 1, 0, 0, 0],
    [1, 5, 1, 0, 1],
    [0, 1, 5, 1, 0],
    [0, 0, 1, 5, 1],
    [0, 1, 0, 1, 5]
]

calculate_bandwidth_and_show_steps(sample_matrix)