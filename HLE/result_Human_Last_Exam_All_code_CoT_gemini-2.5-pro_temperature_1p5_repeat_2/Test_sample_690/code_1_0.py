import math

def find_leftmost_flawed(row, search_range_end):
    """Implements the flawed binary search for the leftmost non-zero element."""
    left = 0
    right = search_range_end
    # Initialize with an index that implies no non-zero was found to the left.
    leftmost_col = -1

    while left <= right:
        mid = (left + right) // 2
        if row[mid] == 0:
            left = mid + 1
        else:  # non-zero
            leftmost_col = mid
            right = mid - 1
    return leftmost_col

def find_rightmost_flawed(row, search_range_start):
    """Implements the flawed binary search for the rightmost non-zero element."""
    n = len(row)
    left = search_range_start
    right = n - 1
    # Initialize with an index that implies no non-zero was found to the right.
    rightmost_col = -1

    while left <= right:
        mid = (left + right) // 2
        if row[mid] == 0:
            right = mid - 1
        else:  # non-zero
            rightmost_col = mid
            left = mid + 1
    return rightmost_col

def flawed_bandwidth_algo(matrix):
    """Implements the flawed algorithm described in the problem."""
    print("--- Running Flawed Algorithm ---")
    n = len(matrix)
    if n == 0:
        return 0
    
    overall_bandwidth = 0
    print("Row | Leftmost Col | Rightmost Col | Row Bandwidth | Equation")
    print("----|--------------|---------------|---------------|---------------------------")

    for i in range(n):
        row = matrix[i]
        
        # Find leftmost non-zero in row[0...i]
        j_left_candidate = find_leftmost_flawed(row, i)
        j_left = j_left_candidate if j_left_candidate != -1 else i

        # Find rightmost non-zero in row[i...n-1]
        j_right_candidate = find_rightmost_flawed(row, i)
        j_right = j_right_candidate if j_right_candidate != -1 else i

        left_dist = i - j_left
        right_dist = j_right - i
        
        row_bandwidth = max(left_dist, right_dist)
        overall_bandwidth = max(overall_bandwidth, row_bandwidth)
        
        equation = f"max({i} - {j_left}, {j_right} - {i}) = max({left_dist}, {right_dist})"
        print(f" {i:2} | {j_left:12} | {j_right:13} | {row_bandwidth:13} | {equation}")
        
    print(f"\nFinal Bandwidth from flawed algorithm: {overall_bandwidth}")
    return overall_bandwidth

def correct_bandwidth_algo(matrix):
    """A simple, correct algorithm to find bandwidth using an O(n^2) scan."""
    print("\n--- Running Correct Algorithm (O(n^2) scan) ---")
    n = len(matrix)
    if n == 0:
        return 0
        
    bandwidth = 0
    for i in range(n):
        for j in range(n):
            if matrix[i][j] != 0:
                dist = abs(i - j)
                if dist > bandwidth:
                    bandwidth = dist
    print(f"Correct Bandwidth: {bandwidth}")
    return bandwidth

# A symmetric matrix where the flawed algorithm fails
counter_example_matrix = [
    [1, 0, 1],
    [0, 1, 0],
    [1, 0, 1]
]

print("Analyzing matrix:")
for r in counter_example_matrix:
    print(r)
print("-" * 30)

flawed_bandwidth_algo(counter_example_matrix)
correct_bandwidth_algo(counter_example_matrix)
