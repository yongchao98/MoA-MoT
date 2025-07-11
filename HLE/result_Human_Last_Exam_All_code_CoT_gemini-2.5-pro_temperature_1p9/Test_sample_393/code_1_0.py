import math

def calculate_tverberg_partitions(partition_parts):
    """
    Calculates the number of Tverberg partitions for a given integer partition.

    Args:
      partition_parts: A list of integers [a_1, a_2, ...] representing the
                       sizes of the parts in the partition of n.

    Returns:
      The number of partitions, which is calculated using the multinomial
      coefficient formula: n! / (a_1! * a_2! * ... * a_r!).
    """
    a = partition_parts
    if not all(isinstance(x, int) and x > 0 for x in a):
        print("Error: All partition parts must be positive integers.")
        return None

    n = sum(a)
    r = len(a)

    # Calculate the multinomial coefficient
    try:
        numerator = math.factorial(n)
        denominator = 1
        for part_size in a:
            denominator *= math.factorial(part_size)
        
        result = numerator // denominator
    except ValueError:
        print("Error: Input contains negative numbers, which is not valid for factorials.")
        return None


    # Output the explanation and the step-by-step calculation
    print(f"For the partition of n = {n} into r = {r} parts: {a}")
    print("The number of Tverberg partitions is calculated as follows:")
    
    # Print the symbolic formula
    denom_symbolic_str = " * ".join([f"{p}!" for p in a])
    print(f"  Formula: {n}! / ({denom_symbolic_str})")
    
    # Print the equation with factorial values
    numerator_val_str = str(numerator)
    denominator_vals_str = " * ".join([str(math.factorial(p)) for p in a])
    print(f"= {numerator_val_str} / ({denominator_vals_str})")
    print(f"= {numerator_val_str} / {denominator}")
    print(f"= {result}")
    
    return result

if __name__ == '__main__':
    # As the problem does not provide a specific partition, an example is used.
    # Let's consider a partition of n=8 into r=4 parts.
    # Example partition: [3, 2, 2, 1]
    # For this partition to be valid, the condition max(a_i) <= d+1 must hold.
    # Here, max(a_i) = 3, so we need 3 <= d+1, which means d >= 2.
    example_partition = [3, 2, 2, 1]
    
    calculate_tverberg_partitions(example_partition)
    # The result for this specific example is 3360.
    # If another partition was chosen, the result would change. For example, for [2,2,2,2] it would be 2520.
