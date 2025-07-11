def get_powerset(s):
    """Computes the powerset of a given set."""
    x = len(s)
    powerset = set()
    for i in range(1 << x):
        subset = frozenset({e for j, e in enumerate(s) if (i >> j) & 1})
        powerset.add(subset)
    return powerset

def get_dual_topology(topology, base_set):
    """
    Calculates the dual topology for a given topology on a finite set.
    In a finite space:
    - All subsets are compact.
    - Saturated sets (intersections of open sets) are just open sets.
    - Thus, the closed sub-basis for the dual is the original topology itself.
    """
    
    # The closed basis is formed by finite unions of sub-basis elements.
    # Since the sub-basis is a topology, this is just the topology itself.
    closed_basis = topology

    # The closed sets are arbitrary intersections of the basis elements.
    # Since any topology on a finite set is Alexandroff (closed under
    # arbitrary intersections), the closed sets are again just the topology itself.
    dual_closed_sets = topology

    # The open sets of the dual topology are the complements of its closed sets.
    dual_open_sets = {base_set - s for s in dual_closed_sets}
    
    # We must ensure the result is a set of frozensets for hashing
    return frozenset(dual_open_sets)

def print_topology(t, name):
    """Helper function to print a topology in a readable format."""
    print(f"{name}: {{ {', '.join(sorted(['{' + ', '.join(map(str, sorted(list(s)))) + '}' for s in t]))} }}")

# Define a base set for our topologies
X = frozenset({1, 2, 3})
EMPTYSET = frozenset()

# Define an initial topology T0 on X.
# This one is not self-complementary.
T0 = frozenset({
    EMPTYSET,
    frozenset({1}),
    frozenset({1, 2}),
    X
})

# --- Main Execution ---
print(f"Let's analyze the iteration of the dual topology operator on a finite set X = {{1, 2, 3}}.\n")

# Start with our initial topology T0
print_topology(T0, "T_0 (Initial)")
print("-" * 20)

# First iteration: T1 = d(T0)
T1 = get_dual_topology(T0, X)
print_topology(T1, "T_1 = dual(T_0)")

# Second iteration: T2 = d(T1)
T2 = get_dual_topology(T1, X)
print_topology(T2, "T_2 = dual(T_1)")
print("-" * 20)

# Check for stabilization
print("Checking for stabilization:")
if T2 == T0:
    print("Result: T_2 is identical to T_0.")
    print("The sequence of distinct topologies is (T_0, T_1).")
    num_distinct = len({T0, T1})
    print(f"The number of distinct topologies generated is {num_distinct}.")
else:
    print("Error in logic: T_2 was not equal to T_0.")

# The maximum possible number is when T_1 is not equal to T_0.
max_num = 2
print(f"\nFor any topology on a finite set, the sequence is T_0, T_1, T_0, ...")
print(f"Therefore, the largest possible number of distinct topologies is {max_num}.")
