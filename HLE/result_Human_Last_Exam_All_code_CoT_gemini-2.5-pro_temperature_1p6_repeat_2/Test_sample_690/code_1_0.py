import math

def calculate_bandwidth(matrix):
    """
    Calculates the bandwidth of a square matrix using the specified algorithm.

    The function implements the following logic:
    1. Iterates through each row of the matrix.
    2. For each row `i`, it splits the row at the diagonal.
    3. It uses binary search to find the leftmost non-zero element in the range [0, i].
    4. It uses binary search to find the rightmost non-zero element in the range [i, n-1].
    5. It calculates the row's bandwidth as the maximum distance of these elements from the diagonal.
    6. The final bandwidth is the maximum of all row bandwidths.
    """
    n = len(matrix)
    if n == 0:
        print("Bandwidth = 0")
        return

    # Helper for finding the leftmost non-zero element in a list
    def find_leftmost_nonzero(arr):
        left, right = 0, len(arr) - 1
        leftmost_idx = -1
        while left <= right:
            mid = (left + right) // 2
            if arr[mid] != 0:
                leftmost_idx = mid
                right = mid - 1
            else:
                left = mid + 1
        return leftmost_idx

    # Helper for finding the rightmost non-zero element in a list
    def find_rightmost_nonzero(arr):
        left, right = 0, len(arr) - 1
        rightmost_idx = -1
        while left <= right:
            mid = (left + right) // 2
            if arr[mid] != 0:
                rightmost_idx = mid
                left = mid + 1
            else:
                right = mid - 1
        return rightmost_idx

    overall_bandwidth = 0
    row_bandwidths = []

    print("Analyzing matrix:")
    for row in matrix:
        print(row)
    print("-" * 20)

    for i in range(n):
        # b. Find the leftmost non-zero element up to the diagonal
        left_half = matrix[i][:i+1]
        leftmost_col_local = find_leftmost_nonzero(left_half)
        # The global column index is the same as the local index
        leftmost_col_global = leftmost_col_local

        # c. Find the rightmost non-zero element from the diagonal onwards
        right_half = matrix[i][i:]
        rightmost_col_local = find_rightmost_nonzero(right_half)
        # The global column index is the local index plus the starting offset i
        rightmost_col_global = rightmost_col_local + i if rightmost_col_local != -1 else -1

        # d. Calculate the distance from the diagonal
        dist_left = 0
        if leftmost_col_global != -1:
            dist_left = i - leftmost_col_global

        dist_right = 0
        if rightmost_col_global != -1:
            dist_right = rightmost_col_global - i

        row_bw = max(dist_left, dist_right)
        row_bandwidths.append(row_bw)
        
        # e. Update overall bandwidth
        if row_bw > overall_bandwidth:
            overall_bandwidth = row_bw
            
        # Print intermediate steps for clarity
        print(f"Row {i}:")
        left_str = f"Leftmost nz at col {leftmost_col_global}" if leftmost_col_global != -1 else "No left nz"
        right_str = f"Rightmost nz at col {rightmost_col_global}" if rightmost_col_global != -1 else "No right nz"
        print(f"  {left_str}, {right_str}")
        print(f"  Distances from diagonal: |{i}-{leftmost_col_global}|={dist_left}, |{rightmost_col_global}-{i}|={dist_right}")
        print(f"  Row bandwidth = {row_bw}")


    print("-" * 20)
    # 6. Output the value of the bandwidth and the 'equation'
    equation_str = "max(" + ", ".join(map(str, row_bandwidths)) + ")"
    print(f"Final Equation: Bandwidth = {equation_str}")
    print(f"Result: {overall_bandwidth}")


if __name__ == '__main__':
    # Example symmetric matrix from Wikipedia
    symmetric_matrix = [
        [3, 1, 0, 0, 0],
        [1, 3, 1, 0, 0],
        [0, 1, 3, 1, 0],
        [0, 0, 1, 3, 1],
        [0, 0, 0, 1, 3]
    ]
    
    # Another example with a larger bandwidth
    symmetric_matrix_2 = [
        [5, 0, 2, 0, 0],
        [0, 6, 0, 3, 0],
        [2, 0, 7, 0, 4],
        [0, 3, 0, 8, 0],
        [0, 0, 4, 0, 9]
    ]
    
    # Use one of the examples
    calculate_bandwidth(symmetric_matrix_2)

<<<A>>>