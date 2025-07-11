import math
from collections import Counter

def count_tverberg_partitions(partition_sizes: list[int]):
    """
    Calculates the number of Tverberg partitions for a set of n points,
    given that all partitions have parts of the specified sizes.

    The problem is equivalent to counting the number of ways to partition a set of
    n elements into subsets of given sizes a_1, a_2, ..., a_r.

    The formula is: n! / (a_1! * a_2! * ... * a_r!) / (c_1! * c_2! * ... * c_k!)
    where c_j is the count of each distinct part size.
    """
    if not partition_sizes:
        print("The list of partition sizes cannot be empty.")
        return

    n = sum(partition_sizes)
    r = len(partition_sizes)

    # Calculate the multinomial coefficient part: n! / (a_1! * a_2! * ... * a_r!)
    # We do this in log space to avoid potential overflow with large numbers, though for
    # this specific calculation, direct computation is often fine.
    # However, direct computation is easier to show in the final string.
    
    numerator = math.factorial(n)
    
    parts_factorial_product = 1
    for part in partition_sizes:
        if part <= 0:
            print("Error: Partition sizes must be positive integers.")
            return
        parts_factorial_product *= math.factorial(part)

    # Count frequencies of each part size to handle indistinguishable groups
    size_counts = Counter(partition_sizes)
    
    counts_factorial_product = 1
    for count in size_counts.values():
        counts_factorial_product *= math.factorial(count)
        
    # Final result
    try:
        if parts_factorial_product * counts_factorial_product == 0:
             raise ZeroDivisionError
        result = numerator // (parts_factorial_product * counts_factorial_product)
    except ZeroDivisionError:
        print("Error: Division by zero. This might be caused by non-positive partition sizes.")
        return

    # --- Formatting the output equation ---
    
    # Numerator string
    num_str = f"{n}!"
    
    # Denominator part 1: (a_1! * a_2! * ...)
    parts_fact_str = " * ".join([f"{p}!" for p in partition_sizes])
    
    # Denominator part 2: (c_1! * c_2! * ...)
    counts_fact_str = " * ".join([f"{c}!" for c in size_counts.values()])

    print("The number of Tverberg partitions is calculated as follows:")
    print(f"Number = {num_str} / ({parts_fact_str}) / ({counts_fact_str})")
    
    num_val = numerator
    parts_fact_val = parts_factorial_product
    counts_fact_val = counts_factorial_product
    
    print(f"Number = {num_val} / ({parts_fact_val}) / ({counts_fact_val})")
    print(f"Number = {result}")


if __name__ == '__main__':
    # Example: A partition of n=8 into parts of sizes 3, 2, 2, 1
    # n = 3+2+2+1 = 8
    # r = 4
    # The list of parts is a = [3, 2, 2, 1]
    # The distinct sizes are 3 (appears 1 time), 2 (appears 2 times), 1 (appears 1 time)
    # The counts are c = [1, 2, 1]
    # Formula: 8! / (3! * 2! * 2! * 1!) / (1! * 2! * 1!)
    # 40320 / (6 * 2 * 2 * 1) / (1 * 2 * 1) = 40320 / 24 / 2 = 1680 / 2 = 840
    
    example_partition = [3, 2, 2, 1]
    count_tverberg_partitions(example_partition)