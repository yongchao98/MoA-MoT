import math
from collections import Counter

def get_partitions_recursive(n, min_part=1):
    """
    Generate all partitions of an integer n as sorted lists.
    """
    if n == 0:
        yield []
        return
    for i in range(min_part, n + 1):
        for p in get_partitions_recursive(n - i, i):
            yield [i] + p

def solve():
    """
    Solves the problem by calculating the number of non-trivial conjugacy classes in A_5.
    
    The rank of H_c^2(Y, Q) for a crepant resolution Y of C^3/G (with G=A_5)
    is equal to the number of non-trivial conjugacy classes of G.
    
    The script calculates this number by:
    1. Finding all partitions of 5, which correspond to conjugacy classes in S_5.
    2. Identifying which partitions correspond to even permutations (classes in A_5).
    3. Checking the splitting condition for each even class to find the number of A_5 classes.
    4. Subtracting 1 for the trivial (identity) class.
    """
    n = 5
    
    partitions = list(get_partitions_recursive(n))

    total_classes_in_A5 = 0
    print("Finding the number of conjugacy classes of A_5:")
    print("-" * 60)
    print(f"{'Partition of 5':<18} {'Permutation Type':<12} {'Splits in A_5?':<18} {'# Classes in A_5'}")
    print("-" * 60)
    
    class_contributions = {}

    for p in partitions:
        # A permutation is even if n - number_of_cycles is even.
        # Here, number_of_cycles is len(p).
        is_even = (n - len(p)) % 2 == 0
        
        perm_type = 'Even' if is_even else 'Odd'
        
        if not is_even:
            print(f"{str(p):<18} {perm_type:<12} {'-':<18} {'0'}")
            continue

        # Check the splitting condition for an S_n class in A_n:
        # The cycle lengths must be distinct and odd.
        counts = Counter(p)
        is_distinct = all(v == 1 for v in counts.values())
        is_all_odd = all(part % 2 != 0 for part in p)
        
        splits = is_distinct and is_all_odd
        
        contribution = 2 if splits else 1
        total_classes_in_A5 += contribution
        class_contributions[str(p)] = contribution
        
        print(f"{str(p):<18} {perm_type:<12} {str(splits):<18} {contribution}")

    print("-" * 60)
    print(f"Total number of conjugacy classes in A_5: {total_classes_in_A5}")
    
    identity_partition_str = str([1, 1, 1, 1, 1])
    num_nontrivial_classes = total_classes_in_A5 - class_contributions[identity_partition_str]
    
    print(f"The identity class corresponds to partition {identity_partition_str}.")
    print(f"Number of non-trivial conjugacy classes = {total_classes_in_A5} - {class_contributions[identity_partition_str]} = {num_nontrivial_classes}")
    print("\nFinal Calculation:")
    
    rank = num_nontrivial_classes
    print(f"rank(H^2_c(Y, Q)) = (Number of conjugacy classes of A_5) - 1")
    print(f"rank(H^2_c(Y, Q)) = {total_classes_in_A5} - 1 = {rank}")
    
    print(f"\nThus, the rank of H^2_c(Y, Q) is {rank}.")
    print(f'<<<{rank}>>>')

solve()