import math

def count_multiset_combinations(n, k):
    """
    Calculates the number of ways to choose k items from a set of n items with replacement.
    This is the multiset coefficient, given by C(n + k - 1, k).
    """
    if n == 0 and k > 0:
        return 0
    return math.comb(n + k - 1, k)

def solve_del_pezzo_problem():
    """
    Calculates the number of isomorphism classes of del Pezzo surfaces of degree 5
    over Q with good reduction everywhere except possibly at the prime 2.
    """
    # Number of number fields unramified outside {2} for degrees 1 to 5.
    # This data is based on number theory results and databases (e.g., LMFDB).
    field_counts = {
        1: 1,  # Q
        2: 3,  # Q(i), Q(sqrt(2)), Q(sqrt(-2))
        3: 0,
        4: 9,
        5: 0,
    }
    
    # We count the number of quintic etale algebras by considering partitions of 5.
    
    # Partition 5: A single field of degree 5.
    # We need to choose 1 field of degree 5 from the available pool.
    n5 = field_counts[5]
    k5 = 1
    count_p5 = count_multiset_combinations(n5, k5)
    
    # Partition 4+1: One field of degree 4 and one field of degree 1.
    n4, n1 = field_counts[4], field_counts[1]
    k4, k1 = 1, 1
    count_p41 = count_multiset_combinations(n4, k4) * count_multiset_combinations(n1, k1)

    # Partition 3+2: One field of degree 3 and one of degree 2.
    n3, n2 = field_counts[3], field_counts[2]
    k3, k2 = 1, 1
    count_p32 = count_multiset_combinations(n3, k3) * count_multiset_combinations(n2, k2)

    # Partition 3+1+1: One field of degree 3 and two of degree 1.
    n3, n1 = field_counts[3], field_counts[1]
    k3, k1 = 1, 2
    count_p311 = count_multiset_combinations(n3, k3) * count_multiset_combinations(n1, k1)

    # Partition 2+2+1: Two fields of degree 2 and one of degree 1.
    n2, n1 = field_counts[2], field_counts[1]
    k2, k1 = 2, 1
    count_p221 = count_multiset_combinations(n2, k2) * count_multiset_combinations(n1, k1)
    
    # Partition 2+1+1+1: One field of degree 2 and three of degree 1.
    n2, n1 = field_counts[2], field_counts[1]
    k2, k1 = 1, 3
    count_p2111 = count_multiset_combinations(n2, k2) * count_multiset_combinations(n1, k1)
    
    # Partition 1+1+1+1+1: Five fields of degree 1.
    n1 = field_counts[1]
    k1 = 5
    count_p11111 = count_multiset_combinations(n1, k1)

    total_count = (count_p5 + count_p41 + count_p32 + count_p311 + 
                   count_p221 + count_p2111 + count_p11111)

    print("The total number of isomorphism classes is the sum of counts for each partition of 5:")
    print(f"Partition 5: {count_p5}")
    print(f"Partition 4+1: {count_p41}")
    print(f"Partition 3+2: {count_p32}")
    print(f"Partition 3+1+1: {count_p311}")
    print(f"Partition 2+2+1: {count_p221}")
    print(f"Partition 2+1+1+1: {count_p2111}")
    print(f"Partition 1+1+1+1+1: {count_p11111}")
    print("-" * 30)
    print(f"Total = {count_p5} + {count_p41} + {count_p32} + {count_p311} + {count_p221} + {count_p2111} + {count_p11111} = {total_count}")

solve_del_pezzo_problem()