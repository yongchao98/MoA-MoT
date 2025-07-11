def calculate_bandwidth(matrix):
    """
    Calculates the bandwidth of a square matrix using the described algorithm.
    It prints the step-by-step process.
    """
    n = len(matrix)
    if n == 0:
        print("Matrix is empty, bandwidth is 0.")
        return 0

    print(f"Matrix size: {n}x{n}\n")
    # 1. Initialize the bandwidth to 0.
    overall_bandwidth = 0

    # 2. Iterate through each row of the matrix.
    for i in range(n):
        print(f"--- Processing Row {i} ---")
        
        # 3.b Find the leftmost non-zero element
        # Binary search between column 0 and column i
        left = 0
        right = i
        leftmost_col = -1
        while left <= right:
            mid = (left + right) // 2
            if matrix[i][mid] != 0:
                leftmost_col = mid
                right = mid - 1
            else:
                left = mid + 1
        
        if leftmost_col != -1:
            print(f"Leftmost non-zero element found at column: {leftmost_col}")
        else:
            print(f"No non-zero element found to the left of or on the diagonal.")

        # 3.c Find the rightmost non-zero element
        # Binary search between column i and column n-1
        left = i
        right = n - 1
        rightmost_col = -1
        while left <= right:
            mid = (left + right) // 2
            if matrix[i][mid] != 0:
                rightmost_col = mid
                left = mid + 1
            else:
                right = mid - 1
        
        if rightmost_col != -1:
            print(f"Rightmost non-zero element found at column: {rightmost_col}")
        else:
            print(f"No non-zero element found to the right of or on the diagonal.")
            
        # 3.d Calculate the distance from the diagonal
        dist_left = 0
        if leftmost_col != -1:
            dist_left = i - leftmost_col
            print(f"Distance from diagonal for leftmost: {i} - {leftmost_col} = {dist_left}")

        dist_right = 0
        if rightmost_col != -1:
            dist_right = rightmost_col - i
            print(f"Distance from diagonal for rightmost: {rightmost_col} - {i} = {dist_right}")

        row_bandwidth = max(dist_left, dist_right)
        print(f"Row bandwidth = max({dist_left}, {dist_right}) = {row_bandwidth}")

        # 3.e Update the overall bandwidth
        old_overall_bandwidth = overall_bandwidth
        overall_bandwidth = max(overall_bandwidth, row_bandwidth)
        print(f"Overall bandwidth = max({old_overall_bandwidth}, {row_bandwidth}) = {overall_bandwidth}\n")

    # 6. Output the value of the bandwidth.
    print(f"Final calculated bandwidth: {overall_bandwidth}")
    return overall_bandwidth

if __name__ == '__main__':
    # A sample 5x5 real-valued, symmetric band matrix.
    # The largest distance is for A[0][2] or A[2][0], where |2-0|=2
    # So the expected bandwidth is 2.
    sample_matrix = [
        [9, -3, 6, 0, 0],
        [-3, 7, 1, 0, 0],
        [6, 1, 5, 2, 0],
        [0, 0, 2, 8, 4],
        [0, 0, 0, 4, 3]
    ]
    calculate_bandwidth(sample_matrix)