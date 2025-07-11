def flawed_leftmost_search(row, start, end):
    """
    Implements the flawed binary search for the leftmost non-zero element
    in the subarray row[start:end+1] as described in the problem.
    """
    left = start
    right = end
    # Initialize to a value indicating not found.
    # If no non-zero is found, this remains, and we can handle it.
    leftmost_col = -1

    while left <= right:
        mid = (left + right) // 2
        if row[mid] != 0:
            leftmost_col = mid
            right = mid - 1
        else:  # row[mid] == 0
            left = mid + 1
    return leftmost_col


def flawed_rightmost_search(row, start, end):
    """
    Implements the flawed binary search for the rightmost non-zero element
    in the subarray row[start:end+1] as described in the problem.
    """
    left = start
    right = end
    rightmost_col = -1

    while left <= right:
        mid = (left + right) // 2
        if row[mid] != 0:
            rightmost_col = mid
            left = mid + 1
        else:  # row[mid] == 0
            right = mid - 1
    return rightmost_col


def compute_bandwidth_with_flawed_algorithm(matrix):
    """
    Implements the full flawed algorithm and prints its step-by-step execution.
    """
    n = len(matrix)
    if n == 0:
        return 0

    bandwidth = 0
    print("Executing the algorithm from the problem description:")
    print("=" * 50)

    for i in range(n):
        print(f"Processing Row {i}: {matrix[i]}")

        # Find leftmost non-zero in row i, from column 0 to i
        leftmost_col = flawed_leftmost_search(matrix[i], 0, i)
        if leftmost_col == -1:
            dist_left = 0
            print(f"  Leftmost non-zero (cols 0-{i}): Not found. Distance = 0")
        else:
            dist_left = i - leftmost_col
            print(f"  Leftmost non-zero (cols 0-{i}): Found at col {leftmost_col}. Distance = |{i} - {leftmost_col}| = {dist_left}")

        # Find rightmost non-zero in row i, from column i to n-1
        rightmost_col = flawed_rightmost_search(matrix[i], i, n - 1)
        if rightmost_col == -1:
            dist_right = 0
            print(f"  Rightmost non-zero (cols {i}-{n-1}): Not found. Distance = 0")
        else:
            dist_right = rightmost_col - i
            print(f"  Rightmost non-zero (cols {i}-{n-1}): Found at col {rightmost_col}. Distance = |{rightmost_col} - {i}| = {dist_right}")

        # Calculate row bandwidth and update overall bandwidth
        row_bandwidth = max(dist_left, dist_right)
        print(f"  Row {i} Bandwidth Equation: max({dist_left}, {dist_right}) = {row_bandwidth}")

        if row_bandwidth > bandwidth:
            print(f"  Updating overall bandwidth from {bandwidth} to {row_bandwidth}")
            bandwidth = row_bandwidth
        print("-" * 20)
    
    return bandwidth


def compute_bandwidth_correctly(matrix):
    """A correct (but less efficient) algorithm for comparison."""
    n = len(matrix)
    if n == 0:
        return 0
    max_dist = 0
    for i in range(n):
        for j in range(n):
            if matrix[i][j] != 0:
                dist = abs(i-j)
                if dist > max_dist:
                    max_dist = dist
    return max_dist


# A symmetric matrix where the flawed algorithm fails
test_matrix = [
    [2, 0, 3],
    [0, 1, 0],
    [3, 0, 5]
]

print("Demonstrating the flaw with a test matrix:")
print(f"Test Matrix = {test_matrix}\n")

flawed_result = compute_bandwidth_with_flawed_algorithm(test_matrix)
correct_result = compute_bandwidth_correctly(test_matrix)

print("\n" + "=" * 50)
print("Final Comparison:")
print(f"Result from the described algorithm: {flawed_result}")
print(f"Correct bandwidth: {correct_result}")
print("Conclusion: The algorithm is flawed and produced an incorrect result.")
