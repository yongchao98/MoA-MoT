import math
from collections import Counter

def calculate_tverberg_partitions():
    """
    Calculates the number of Tverberg partitions for a given integer partition.

    The problem reduces to a combinatorial task of counting the number of ways to
    partition a set of n elements into unlabeled subsets of specified sizes.
    """
    # The problem is stated for any valid partition a_1, ..., a_r.
    # Since no specific partition is provided, we will use an example
    # to demonstrate the calculation.
    # Let the partition of n be given by the parts in list 'a'.
    a = [3, 2, 2, 1]
    
    # 1. Calculate n, the total number of points.
    n = sum(a)
    
    # 2. Calculate the numerator: n!
    n_factorial = math.factorial(n)
    
    # 3. Calculate the first part of the denominator: product of a_i!
    # This accounts for the number of ways to arrange elements within each part.
    part_factorials = [math.factorial(part_size) for part_size in a]
    part_factorials_prod = math.prod(part_factorials)
    
    # 4. Calculate the second part of the denominator: product of c_s!
    # This accounts for the overcounting of partitions with parts of the same size.
    # First, count the frequency of each part size.
    size_counts = Counter(a)
    # Then, calculate the product of the factorials of these counts.
    count_factorials = [math.factorial(count) for count in size_counts.values()]
    count_factorials_prod = math.prod(count_factorials)
    
    # 5. Calculate the final number of partitions.
    # Using integer division as the result must be an integer.
    num_partitions = n_factorial // (part_factorials_prod * count_factorials_prod)

    # 6. Print the details of the calculation as requested.
    print(f"Given the partition a = {a}")
    print(f"The total number of points is n = {' + '.join(map(str, a))} = {n}")
    print("\nThe number of Tverberg's partitions is given by the formula:")
    print("N = n! / ( (a_1! * a_2! * ...) * (c_1! * c_2! * ...) )")
    print("where a_i are the sizes of the parts and c_s are the frequencies of those sizes.")
    
    part_fact_str = " * ".join([f"{x}!" for x in a])
    count_list_str = " * ".join([f"{c}!" for c in size_counts.values()])
    
    print("\nCalculation steps:")
    print(f"Numerator: n! = {n}! = {n_factorial}")
    
    print("\nDenominator Part 1 (Product of factorials of part sizes):")
    print(f"  {part_fact_str} = {' * '.join(map(str, part_factorials))} = {part_factorials_prod}")

    print("\nDenominator Part 2 (Product of factorials of size frequencies):")
    print(f"  The frequencies of part sizes {list(size_counts.keys())} are {list(size_counts.values())}.")
    print(f"  {count_list_str} = {' * '.join(map(str, count_factorials))} = {count_factorials_prod}")
    
    print("\nFinal Result:")
    print(f"Number of partitions = {n_factorial} / ({part_factorials_prod} * {count_factorials_prod})")
    print(f"                     = {n_factorial} / {part_factorials_prod * count_factorials_prod}")
    print(f"                     = {num_partitions}")


# Execute the function
calculate_tverberg_partitions()