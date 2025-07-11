import math

def partitions(n, max_val=None):
    """
    Generates all integer partitions of n.
    An integer partition is a way of writing n as a sum of positive integers.
    The generator yields partitions as tuples in reverse lexicographical order.
    For example, for n=4, it yields (4,), (3, 1), (2, 2), (2, 1, 1), (1, 1, 1, 1).
    """
    if max_val is None:
        max_val = n
    if n == 0:
        yield ()
        return
    # 'i' is the first and largest part of the partition.
    for i in range(min(n, max_val), 0, -1):
        # The rest of the partition must sum to n-i, and its parts
        # must be no larger than i to maintain descending order.
        for p in partitions(n - i, i):
            yield (i,) + p

def get_irrep_dimension(n, partition):
    """
    Calculates the dimension of an irreducible representation of S_n
    corresponding to the given partition, using the Hook Length Formula.
    """
    # The numerator in the formula is n!
    numerator = math.factorial(n)
    
    # The denominator is the product of all hook lengths.
    denominator = 1
    
    # We iterate through each cell of the Young diagram for the partition.
    # The row index is 'i' and the column index is 'j'.
    for i, row_len in enumerate(partition):
        for j in range(row_len):
            # Calculate the hook length for the cell at (i, j).
            
            # Count cells to the right in the same row, plus the cell itself.
            hooks_in_row = row_len - j
            
            # Count cells below in the same column, plus the cell itself.
            hooks_in_col = 0
            for k in range(i, len(partition)):
                # A cell exists at row k, column j if the length of row k is at least j+1.
                if partition[k] > j:
                    hooks_in_col += 1
                else:
                    # Since partition parts are sorted in descending order,
                    # no subsequent row will have a cell in this column.
                    break
            
            # The hook length is the sum of hooks in the row and column,
            # minus 1 because the cell itself is counted in both.
            hook_length = hooks_in_row + hooks_in_col - 1
            denominator *= hook_length
    
    # The dimension is guaranteed to be an integer. Python's arbitrary-precision
    # integers handle the large numbers, and integer division gives the exact result.
    return numerator // denominator

def main():
    """
    Main function to solve the problem.
    """
    n = 25
    dimension_threshold = 500000
    
    # This will count the number of irreps with dimension below the threshold.
    count = 0
    
    # Generate all partitions of n.
    partition_generator = partitions(n)
    
    for p in partition_generator:
        # For each partition, calculate the dimension of the corresponding irrep.
        dimension = get_irrep_dimension(n, p)
        
        # Check if the dimension is strictly less than the threshold.
        if dimension < dimension_threshold:
            count += 1
            
    # Print the final result, including the numbers from the problem statement.
    print(f"The number of irreducible representations of the symmetric group S_{n} "
          f"that have a dimension strictly less than {dimension_threshold} is {count}.")

if __name__ == "__main__":
    main()
<<<106>>>