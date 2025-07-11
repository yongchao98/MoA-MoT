import math

def solve():
    """
    Calculates the number of irreducible representations of the symmetric group S_25
    that have a dimension strictly less than 500,000.
    """
    n = 25
    limit = 500_000

    # Pre-calculate n! as it's a constant for all partitions.
    n_factorial = math.factorial(n)

    memo_partitions = {}
    def generate_partitions(num, max_val=None):
        """
        A generator for all integer partitions of 'num' in descending order.
        For example, for n=4, it yields [4], [3,1], [2,2], [2,1,1], [1,1,1,1].
        """
        if max_val is None:
            max_val = num

        # Base case for the recursion.
        if num == 0:
            yield []
            return

        # Recursively find partitions.
        for i in range(min(num, max_val), 0, -1):
            for p_rest in generate_partitions(num - i, i):
                yield [i] + p_rest

    memo_dimension = {}
    def calculate_dimension(partition):
        """
        Calculates the dimension of an irrep for a given partition using the hook-length formula.
        Memoization is used to cache results for previously seen partitions (though for n=25, all partitions are unique).
        """
        # Partitions are lists, which are not hashable. Convert to tuple for memoization key.
        partition_tuple = tuple(partition)
        if partition_tuple in memo_dimension:
            return memo_dimension[partition_tuple]

        product_of_hooks = 1
        # Iterate over each cell (r, c) in the Young diagram (0-indexed).
        for r, row_length in enumerate(partition):
            for c in range(row_length):
                # Number of cells to the right in the same row.
                hooks_right = row_length - (c + 1)
                
                # Number of cells below in the same column.
                hooks_down = 0
                for r_below in range(r + 1, len(partition)):
                    if partition[r_below] > c:
                        hooks_down += 1
                
                # Hook length is the cell itself (1) + right hooks + down hooks.
                hook_length = 1 + hooks_right + hooks_down
                product_of_hooks *= hook_length
        
        # Dimension = n! / product of all hook lengths. Use integer division.
        dimension = n_factorial // product_of_hooks
        
        memo_dimension[partition_tuple] = dimension
        return dimension

    count = 0
    # Iterate through all partitions of n.
    for p in generate_partitions(n):
        dim = calculate_dimension(p)
        if dim < limit:
            count += 1
    
    # Format the final output as a sum of 1s, as requested.
    if count == 0:
        print("0")
    else:
        # Create a string like "1 + 1 + ... + 1 = count"
        sum_of_ones_str = " + ".join(["1"] * count)
        print(f"{sum_of_ones_str} = {count}")

solve()
<<<70>>>