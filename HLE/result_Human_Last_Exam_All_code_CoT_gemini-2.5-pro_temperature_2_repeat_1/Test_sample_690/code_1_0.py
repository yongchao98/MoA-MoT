import math

def find_leftmost_nonzero(row, search_range):
    """
    Performs binary search to find the leftmost non-zero element in a given range of a row.
    Returns column index (0-based) or -1 if not found.
    """
    left, right = search_range
    leftmost_col = -1
    initial_left, initial_right = left, right
    
    while left <= right:
        mid = (left + right) // 2
        if row[mid] != 0:
            leftmost_col = mid
            right = mid - 1
        else:
            left = mid + 1
            
    print(f"    - Leftmost non-zero search in columns {initial_left+1}-{initial_right+1}: ", end="")
    if leftmost_col != -1:
        print(f"found at column {leftmost_col + 1}.")
    else:
        print(f"none found.")
    return leftmost_col

def find_rightmost_nonzero(row, search_range):
    """
    Performs binary search to find the rightmost non-zero element in a given range of a row.
    Returns column index (0-based) or -1 if not found.
    """
    left, right = search_range
    rightmost_col = -1
    initial_left, initial_right = left, right
    
    while left <= right:
        mid = (left + right) // 2
        if row[mid] != 0:
            rightmost_col = mid
            left = mid + 1
        else:
            right = mid - 1
            
    print(f"    - Rightmost non-zero search in columns {initial_left+1}-{initial_right+1}: ", end="")
    if rightmost_col != -1:
        print(f"found at column {rightmost_col + 1}.")
    else:
        print(f"none found.")
    return rightmost_col

def compute_bandwidth_from_description(matrix):
    """
    Implements the algorithm described in the problem to compute matrix bandwidth.
    """
    n = len(matrix)
    if n == 0:
        return 0
        
    overall_bandwidth = 0
    print(f"--- Running Algorithm on a {n}x{n} Matrix ---")
    print(f"Initial overall bandwidth = {overall_bandwidth}\n")

    for i in range(n):
        row_index_1_based = i + 1
        print(f"Processing row {row_index_1_based}:")
        
        row = matrix[i]
        
        leftmost_col = find_leftmost_nonzero(row, (0, i))
        rightmost_col = find_rightmost_nonzero(row, (i, n - 1))
        
        dist_left = 0
        if leftmost_col != -1:
            dist_left = i - leftmost_col
            print(f"    - Left distance = |row {row_index_1_based} - col {leftmost_col + 1}| = |{row_index_1_based} - {leftmost_col + 1}| = {dist_left}")
        else:
            print(f"    - Left distance = 0")

        dist_right = 0
        if rightmost_col != -1:
            dist_right = rightmost_col - i
            print(f"    - Right distance = |col {rightmost_col + 1} - row {row_index_1_based}| = |{rightmost_col + 1} - {row_index_1_based}| = {dist_right}")
        else:
            print(f"    - Right distance = 0")

        row_bandwidth = max(dist_left, dist_right)
        print(f"    - Row bandwidth = max({dist_left}, {dist_right}) = {row_bandwidth}")
        
        if row_bandwidth > overall_bandwidth:
            print(f"    - New maximum bandwidth found.")
            overall_bandwidth = row_bandwidth
        
        print(f"  Current overall bandwidth = {overall_bandwidth}\n")

    print(f"--- Final Result ---")
    print(f"Final calculated bandwidth: {overall_bandwidth}")
    return overall_bandwidth

# A symmetric 4x4 matrix whose largest distance from diagonal is 2.
# e.g., for A[1,3]=31 (0-indexed: matrix[0][2]), the distance is |0-2| = 2
symmetric_matrix_with_k_2 = [
    [11, 21, 31,  0],
    [21, 22, 32, 42],
    [31, 32, 33, 43],
    [ 0, 42, 43, 44],
] 

compute_bandwidth_from_description(symmetric_matrix_with_k_2)
