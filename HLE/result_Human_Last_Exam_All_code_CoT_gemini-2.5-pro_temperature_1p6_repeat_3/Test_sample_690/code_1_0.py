import numpy as np

def find_leftmost_bs(row, left_idx, right_idx):
    """
    Implements the flawed binary search from the description to find the
    leftmost non-zero element in a given 0-indexed row segment.
    - row: the entire row (list or numpy array)
    - left_idx: the 0-indexed start of the search range
    - right_idx: the 0-indexed end of the search range
    Returns the 0-indexed column of the found element, or -1 if not found.
    """
    ans_j = -1
    l, r = left_idx, right_idx
    
    while l <= r:
        mid = (l + r) // 2
        if row[mid] == 0:
            l = mid + 1
        else:
            ans_j = mid
            r = mid - 1
            
    return ans_j

def find_rightmost_bs(row, left_idx, right_idx):
    """
    Implements the flawed binary search from the description to find the
    rightmost non-zero element in a given 0-indexed row segment.
    - row: the entire row (list or numpy array)
    - left_idx: the 0-indexed start of the search range
    - right_idx: the 0-indexed end of the search range
    Returns the 0-indexed column of the found element, or -1 if not found.
    """
    ans_j = -1
    l, r = left_idx, right_idx
    
    while l <= r:
        mid = (l + r) // 2
        if row[mid] == 0:
            r = mid - 1
        else:
            ans_j = mid
            l = mid + 1
            
    return ans_j

def calculate_bandwidth_flawed(matrix):
    """Implements the flawed algorithm from the problem description."""
    n = len(matrix)
    if n == 0:
        return 0
    
    overall_bandwidth = 0
    for i in range(n): # iterate through each row index i
        row = matrix[i]
        
        # Search for leftmost non-zero from column 0 to i
        j_left = find_leftmost_bs(row, 0, i)
        
        # Search for rightmost non-zero from column i to n-1
        j_right = find_rightmost_bs(row, i, n-1)
        
        dist_left = 0
        if j_left != -1:
            dist_left = i - j_left
            
        dist_right = 0
        if j_right != -1:
            dist_right = j_right - i
        
        row_bandwidth = max(dist_left, dist_right)
        
        if row_bandwidth > overall_bandwidth:
            overall_bandwidth = row_bandwidth
            
    return overall_bandwidth

def calculate_bandwidth_correct(matrix):
    """A correct brute-force O(n^2) algorithm for comparison."""
    n = len(matrix)
    if n == 0:
        return 0
        
    max_dist = 0
    for i in range(n):
        for j in range(n):
            if matrix[i][j] != 0:
                dist = abs(i - j)
                if dist > max_dist:
                    max_dist = dist
    return max_dist

# --- Demonstration ---
# Create a symmetric matrix where the flawed algorithm will fail.
# Non-zeros are at (5,5), (5,7), (7,5) using 1-based indexing
# Or (4,4), (4,6), (6,4) using 0-based indexing
n = 8
# The matrix has non-zeros at M[4][4], M[4][6], M[6][4].
# For row i=4, the non-zeros are at columns 4 and 6. This is not monotonic.
M = np.zeros((n, n))
M[4, 4] = 1.0 
M[4, 6] = 1.0
M[6, 4] = 1.0

# The true bandwidth is max(|4-4|, |4-6|, |6-4|) = 2.
correct_bw = calculate_bandwidth_correct(M)

# The flawed algorithm will fail to find the non-zero at (4,6) because
# its binary search for the rightmost element in row 4 starts at (4,4), which
# is non-zero, but the search doesn't proceed further right to find the true rightmost.
flawed_bw = calculate_bandwidth_flawed(M)

print("Demonstration of the algorithm's flaw on a symmetric matrix:")
print(f"Matrix (non-zero elements):\n{np.argwhere(M)}")
print(f"\nCorrect Bandwidth (from O(n^2) algorithm): {correct_bw}")
print(f"Algorithm's Output: {flawed_bw}")

if correct_bw != flawed_bw:
    print("\nThe algorithm produced an incorrect result.")
else:
    print("\nThe algorithm produced a correct result for this specific case (by chance).")
