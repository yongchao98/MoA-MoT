import math

# Step 1: Set up the constants for the problem.
n = 25
limit = 500_000

# Using a dictionary for memoization to speed up partition generation.
memo_partitions = {}

def generate_partitions(num, max_part=None):
    """
    Generates all non-increasing integer partitions of a number 'num'.
    Partitions are returned as a list of lists.
    Uses memoization to avoid recomputing partitions for the same number.
    """
    if max_part is None:
        max_part = num
    
    # Check if the result is already memoized.
    if (num, max_part) in memo_partitions:
        return memo_partitions[(num, max_part)]

    # Base case: The only partition of 0 is an empty list.
    if num == 0:
        return [[]]
    
    partitions_list = []
    # Iterate from min(num, max_part) down to 1 to find parts of the partition.
    for part_val in range(min(num, max_part), 0, -1):
        # Recursively find partitions for the remainder.
        for sub_partition in generate_partitions(num - part_val, part_val):
            partitions_list.append([part_val] + sub_partition)
    
    # Memoize the result before returning.
    memo_partitions[(num, max_part)] = partitions_list
    return partitions_list

def calculate_dimension(partition, n_factorial):
    """
    Calculates the dimension of the irreducible representation for a given partition
    using the Hook Length Formula.
    """
    # Calculate the product of all hook lengths.
    product_of_hooks = 1
    for i, part_len in enumerate(partition):  # i is the row index (0-based)
        for j in range(part_len):  # j is the column index (0-based)
            # Cells to the right in the same row.
            cells_right = part_len - 1 - j
            
            # Cells below in the same column.
            cells_below = 0
            for k in range(i + 1, len(partition)):
                if partition[k] > j:
                    cells_below += 1
                else:
                    # Since partition parts are sorted, we can stop early.
                    break
            
            hook_length = cells_right + cells_below + 1
            product_of_hooks *= hook_length
            
    # The dimension is n! divided by the product of hook lengths.
    # We use integer division `//` as the result is always an integer.
    return n_factorial // product_of_hooks

# Step 2: Pre-calculate n!
n_factorial = math.factorial(n)

# Step 3: Generate all partitions of n.
all_partitions = generate_partitions(n)

# Step 4: Iterate through partitions, calculate dimensions, and count.
count = 0
for p in all_partitions:
    dimension = calculate_dimension(p, n_factorial)
    if dimension < limit:
        count += 1

# Step 5: Print the final result.
print(f"The number of irreducible representations of S_25 with dimension strictly less than 500,000 is {count}.")
<<<44>>>