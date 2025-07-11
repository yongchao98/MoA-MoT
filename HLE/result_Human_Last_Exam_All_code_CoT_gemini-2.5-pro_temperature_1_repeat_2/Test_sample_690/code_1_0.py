def solve():
    """
    Analyzes and implements the bandwidth calculation algorithm.
    """
    # Example symmetric matrix from the problem description analysis.
    # Using 0-based indexing for implementation.
    matrix = [
        [9, 2, 3, 0],
        [2, 9, 0, 4],
        [3, 0, 9, 1],
        [0, 4, 1, 9]
    ]
    n = len(matrix)
    overall_bandwidth = 0
    row_bandwidths = []

    print(f"Analyzing {n}x{n} matrix:")
    for row in matrix:
        print(row)
    print("-" * 20)

    # 1. Initialize bandwidth
    # Using 0-based indexing, i goes from 0 to n-1
    for i in range(n):
        # 3a. Set indices
        # Left search is in range [0, i], right search is in [i, n-1]

        # 3b. Find leftmost non-zero up to diagonal
        left_search_lo, left_search_hi = 0, i
        leftmost_col = -1
        while left_search_lo <= left_search_hi:
            mid = (left_search_lo + left_search_hi) // 2
            if matrix[i][mid] != 0:
                leftmost_col = mid
                left_search_hi = mid - 1
            else:
                left_search_lo = mid + 1

        # 3c. Find rightmost non-zero from diagonal
        right_search_lo, right_search_hi = i, n - 1
        rightmost_col = -1
        while right_search_lo <= right_search_hi:
            mid = (right_search_lo + right_search_hi) // 2
            if matrix[i][mid] != 0:
                rightmost_col = mid
                right_search_lo = mid + 1
            else:
                right_search_hi = mid - 1

        # 3d. Calculate distance from diagonal
        dist_left = i - leftmost_col if leftmost_col != -1 else 0
        dist_right = rightmost_col - i if rightmost_col != -1 else 0
        row_bandwidth = max(dist_left, dist_right)
        row_bandwidths.append(row_bandwidth)
        
        print(f"Row {i+1}: Leftmost non-zero dist = {dist_left}, Rightmost non-zero dist = {dist_right}. Row Bandwidth = {row_bandwidth}")

        # 3e. Update overall bandwidth
        if row_bandwidth > overall_bandwidth:
            overall_bandwidth = row_bandwidth

    # 6. Output the final bandwidth calculation
    print("-" * 20)
    print(f"Final Bandwidth = max({', '.join(map(str, row_bandwidths))}) = {overall_bandwidth}")

solve()