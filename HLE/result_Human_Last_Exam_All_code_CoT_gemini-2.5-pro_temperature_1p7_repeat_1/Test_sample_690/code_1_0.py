import math

def algorithm_from_problem(matrix):
    """
    Implements the flawed algorithm described in the problem statement.
    """
    n = len(matrix)
    if n == 0:
        return 0
    
    # Helper for leftmost non-zero search
    def find_leftmost_in_range(row, start_idx, end_idx):
        left, right = start_idx, end_idx
        leftmost_col = -1
        while left <= right:
            mid = (left + right) // 2
            if row[mid] != 0:
                leftmost_col = mid
                right = mid - 1
            else:
                left = mid + 1
        return leftmost_col

    # Helper for rightmost non-zero search
    def find_rightmost_in_range(row, start_idx, end_idx):
        left, right = start_idx, end_idx
        rightmost_col = -1
        while left <= right:
            mid = (left + right) // 2
            if row[mid] != 0:
                rightmost_col = mid
                left = mid + 1
            else:
                right = mid - 1
        return rightmost_col
    
    overall_bandwidth = 0
    for i in range(n):
        row = matrix[i]
        
        # b. Find the leftmost non-zero element up to the diagonal
        leftmost_j = find_leftmost_in_range(row, 0, i)
        if leftmost_j == -1: # Row is all zeros up to diagonal
            leftmost_j = i
            
        # c. Find the rightmost non-zero element from the diagonal onwards
        rightmost_j = find_rightmost_in_range(row, i, n - 1)
        if rightmost_j == -1: # Row is all zeros from diagonal
             rightmost_j = i
             
        # d. Calculate the distance from the diagonal
        dist_left = i - leftmost_j
        dist_right = rightmost_j - i
        row_bandwidth = max(dist_left, dist_right)
        
        # e. Update overall bandwidth
        if row_bandwidth > overall_bandwidth:
            overall_bandwidth = row_bandwidth
            
    return overall_bandwidth

def correct_bandwidth_algorithm(matrix):
    """
    A simple, correct O(n^2) algorithm to find the bandwidth for comparison.
    """
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

# Create a symmetric matrix where the flawed algorithm will fail
# because of a zero inside the band of non-zeroes.
counter_example_matrix = [
    [5, 0, 2, 0],
    [0, 6, 0, 3],
    [2, 0, 7, 0],
    [0, 3, 0, 8]
]

# Run both algorithms on the matrix
flawed_result = algorithm_from_problem(counter_example_matrix)
correct_result = correct_bandwidth_algorithm(counter_example_matrix)

print("Demonstrating the flaw in the described algorithm:\n")
print("Test Matrix (Symmetric):")
for row in counter_example_matrix:
    print(row)
print("-" * 20)
print(f"Result from the Correct Algorithm (brute-force scan):")
print(f"Bandwidth = {correct_result}")
print(f"Explanation: The furthest non-zero element is A[0,2] or A[1,3], where the distance from the diagonal is abs(0-2)=2 or abs(1-3)=2.\n")

print(f"Result from the Algorithm in the Problem Description:")
print(f"Bandwidth = {flawed_result}")
print("Explanation: The algorithm fails, for example, on row 0 ([5, 0, 2, 0]).")
print("When finding the rightmost non-zero element from the diagonal onwards, the binary search sees the 0 at index 1 and incorrectly concludes the search, missing the 2 at index 2.")
print("This leads to an incorrect final bandwidth calculation.")
print("\nConclusion: The described algorithm is not correct for all symmetric band matrices.")
