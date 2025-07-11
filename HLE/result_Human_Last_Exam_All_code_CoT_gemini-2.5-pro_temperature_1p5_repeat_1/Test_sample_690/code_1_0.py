import math

def find_leftmost_nonzero(row, left, right):
    """
    Performs a binary search to find the index of the leftmost non-zero element
    in the slice row[left:right+1].
    Returns -1 if all elements are zero.
    """
    ans = -1
    while left <= right:
        mid = left + (right - left) // 2
        if row[mid] != 0:
            ans = mid
            right = mid - 1
        else:
            left = mid + 1
    return ans

def find_rightmost_nonzero(row, left, right):
    """
    Performs a binary search to find the index of the rightmost non-zero element
    in the slice row[left:right+1].
    Returns -1 if all elements are zero.
    """
    ans = -1
    while left <= right:
        mid = left + (right - left) // 2
        if row[mid] != 0:
            ans = mid
            left = mid + 1
        else:
            right = mid - 1
    return ans

def compute_bandwidth(matrix):
    """
    Computes the bandwidth of a square matrix using the described algorithm.
    """
    if not matrix or not matrix[0]:
        return 0
    n = len(matrix)
    overall_bandwidth = 0

    print("Analyzing matrix:")
    for row_data in matrix:
        print(row_data)
    print("-" * 20)
    
    # 1. Initialize bandwidth
    max_dist_info = {
        "bandwidth": 0,
        "row": -1,
        "left_dist": 0,
        "right_dist": 0
    }

    # 2. Iterate through each row
    for i in range(n):
        # 3a. Set indices (implicitly done by search function calls)
        current_row = matrix[i]

        # 3b. Find the leftmost non-zero element up to the diagonal
        left_j = find_leftmost_nonzero(current_row, 0, i)
        
        # 3c. Find the rightmost non-zero element from the diagonal onward
        right_j = find_rightmost_nonzero(current_row, i, n - 1)
        
        # 3d. Calculate the distance from the diagonal
        dist_left = i - left_j if left_j != -1 else 0
        dist_right = right_j - i if right_j != -1 else 0
        
        row_bandwidth = max(dist_left, dist_right)
        
        # 3e. Update overall bandwidth if necessary
        if row_bandwidth > max_dist_info["bandwidth"]:
            max_dist_info["bandwidth"] = row_bandwidth
            max_dist_info["row"] = i
            max_dist_info["left_dist"] = dist_left
            max_dist_info["right_dist"] = dist_right
            
    # 6. Output the value of the bandwidth
    # Fulfilling the "output each number in the final equation" request
    # by showing the final calculation that determined the bandwidth.
    bw = max_dist_info['bandwidth']
    r = max_dist_info['row']
    ld = max_dist_info['left_dist']
    rd = max_dist_info['right_dist']
    
    print(f"Maximum distance (bandwidth) found in row {r+1} (0-indexed: {r}).")
    print(f"Calculation for row {r+1}: max(dist_left, dist_right) = max({ld}, {rd}) = {bw}")
    print(f"\nFinal computed bandwidth: {bw}")


# Example usage with a real-valued, symmetric matrix
sample_matrix = [
    [9, 5, 2, 0, 0],
    [5, 8, 4, 0, 0],
    [2, 4, 7, 3, 6],
    [0, 0, 3, 6, 0],
    [0, 0, 6, 0, 5]
]

compute_bandwidth(sample_matrix)

<<<A>>>