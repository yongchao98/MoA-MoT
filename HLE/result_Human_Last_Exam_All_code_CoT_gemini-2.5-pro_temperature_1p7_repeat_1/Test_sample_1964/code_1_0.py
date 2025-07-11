def partition_integers(n_partitions):
    """
    Partitions the set of non-negative integers into n_partitions disjoint "infinite" subsets.
    In this finite simulation, "infinite" just means a list of numbers.
    Each integer k is assigned to the partition k % n_partitions.
    """
    partitions = [[] for _ in range(n_partitions)]
    # Simulate partitioning an infinite set by taking the first 100 integers
    for i in range(1, 101):
        partitions[i % n_partitions].append(i)
    return partitions

def construct_and_verify():
    """
    Illustrates the construction of a Delta-system from a partitioned set.
    """
    # 1. Define a finite root `r`, analogous to a finite set of ordinals.
    r = {0}
    print(f"The chosen finite root is r = {r}\n")

    # 2. Partition the integers (our stand-in for omega \ r).
    # We will create 5 such sets, as an analogue to the omega_1 sets.
    num_sets = 5
    # The partitioned sets P_alpha are disjoint from each other and from r.
    partitions = partition_integers(num_sets)

    # 3. Construct the sets `a_i`, analogous to `a_alpha`.
    # a_i = r U P_i
    A = []
    print("Constructing the collection of sets A = {a_0, a_1, ...}:")
    for i in range(num_sets):
        p_i = set(partitions[i])
        a_i = r.union(p_i)
        A.append(a_i)
        print(f"  a_{i} = {sorted(list(a_i))}")
    print("\nThis collection A is a Delta-system with root r.")

    # 4. Verify the Delta-system property for a pair of sets.
    # The intersection of any two distinct sets a_i, a_j must be the root r.
    i, j = 1, 3
    a_i = A[i]
    a_j = A[j]
    intersection = a_i.intersection(a_j)

    print("\nVerification:")
    print(f"Let's check the intersection of a_{i} and a_{j}:")
    print(f"  a_{i} = {sorted(list(a_i))}")
    print(f"  a_{j} = {sorted(list(a_j))}")
    print(f"  Intersection a_{i} ∩ a_{j} = {sorted(list(intersection))}")

    # Outputting the 'final equation' as requested by the prompt format
    print("\nThe equation for the intersection is:")
    final_eq_str = f"a_{i} ∩ a_{j} = r"
    print(f"  {final_eq_str}")
    
    # We can represent the values in the equation string
    # where the result is the root r.
    result_set = sorted(list(r))
    print(f"Which corresponds to: {sorted(list(a_i))} ∩ {sorted(list(a_j))} = {result_set}")


construct_and_verify()
