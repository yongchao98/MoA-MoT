import math
from collections import Counter

def count_tverberg_partitions(partition_sizes):
    """
    Calculates the number of ways to partition a set of n points into subsets
    with specified sizes.

    This problem is equivalent to counting the number of Tverberg's partitions for
    a specially constructed set of points X, as described in the problem statement.

    Args:
        partition_sizes (list of int): A list of integers [a1, a2, ..., ar]
                                       representing the sizes of the parts in the partition.

    Returns:
        int: The total number of such partitions. Returns -1 if input is invalid.
    """
    if not isinstance(partition_sizes, list) or not all(isinstance(x, int) and x > 0 for x in partition_sizes):
        print("Error: Please provide a list of positive integers for the partition sizes.")
        return -1
    
    if not partition_sizes:
        print("The partition is empty. There is 1 partition of an empty set.")
        return 1

    a = partition_sizes
    r = len(a)
    n = sum(a)

    print(f"Given the partition of n into parts: a = {a}")
    print(f"The total number of points is n = sum(a) = {n}")
    print(f"The number of parts is r = {r}")
    print("-" * 30)

    # The formula for the number of ways to partition a set of n elements into
    # r unlabeled subsets of sizes a_1, ..., a_r is:
    # N = (1 / (c_1! * c_2! * ...)) * (n! / (a_1! * a_2! * ...))
    # where c_i is the number of parts of a certain size.

    # 1. Calculate the multinomial coefficient for ordered partitions.
    # Numerator is n!
    numerator = math.factorial(n)
    
    # Denominator is a_1! * a_2! * ... * a_r!
    denom_a_fact = 1
    for val in a:
        denom_a_fact *= math.factorial(val)
    
    multinomial_coeff = numerator // denom_a_fact

    # 2. Count frequencies of each part size to handle unordered partitions.
    size_counts = Counter(a)
    
    # 3. Calculate the permutation adjustment factor: c_1! * c_2! * ...
    denom_counts_fact = 1
    for count in size_counts.values():
        denom_counts_fact *= math.factorial(count)

    # 4. Calculate the final result
    final_result = multinomial_coeff // denom_counts_fact

    # --- Detailed Output ---
    print("The number of partitions is calculated as follows:")
    print("N = (n! / (a_1! * a_2! * ... * a_r!)) / (c_1! * c_2! * ... * c_k!)")
    
    multinomial_denom_str = " * ".join([f"{val}!" for val in a])
    counts_denom_str = " * ".join([f"{count}!" for count in size_counts.values()])
    print(f"Plugging in the numbers:")
    print(f"N = ({n}! / ({multinomial_denom_str})) / ({counts_denom_str})")
    
    denom_a_fact_vals_str = " * ".join([str(math.factorial(val)) for val in a])
    denom_counts_fact_vals_str = " * ".join([str(math.factorial(count)) for count in size_counts.values()])
    print(f"Calculating the factorial values:")
    print(f"N = ({numerator} / ({denom_a_fact_vals_str})) / ({denom_counts_fact_vals_str})")
    
    print(f"Evaluating the products:")
    print(f"N = ({numerator} / {denom_a_fact}) / {denom_counts_fact}")
    
    print(f"Calculating the multinomial coefficient (for labeled parts):")
    print(f"N = {multinomial_coeff} / {denom_counts_fact}")
    
    print(f"Adjusting for unlabeled parts:")
    print(f"N = {final_result}")
    print("-" * 30)

    return final_result

if __name__ == '__main__':
    # You can replace this list with any partition of n.
    # Example partition a_1, a_2, ..., a_r for a number n.
    example_partition = [3, 3, 1]
    
    final_answer = count_tverberg_partitions(example_partition)
    
    print(f"The total number of Tverberg's partitions for the given structure is: {final_answer}")
    
    # If you need to run this for another partition, just change the list above.
    # For example:
    # count_tverberg_partitions([2, 2])
    # count_tverberg_partitions([4, 1, 1, 1])
