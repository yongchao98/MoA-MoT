import itertools

print("Goal: Find the smallest n for which an n-point space can be reducible (not irreducible).\n")
print("A space is reducible if it's a union of its proper closed subsets.")
print("We will show that n=0 and n=1 spaces are always irreducible, and then construct a reducible space for n=2.\n")

print("--- Analysis for n=2 ---")
# Define a 2-point space X
X = frozenset({0, 1})
print(f"Let the space be X = {set(X)}.")

# Define the discrete topology on X, where every subset is open.
# We use frozensets because they can be elements of a set.
open_sets = {frozenset(), frozenset({0}), frozenset({1}), frozenset({0, 1})}
print("Consider the discrete topology, where every subset is open.")
print(f"Open sets: {[set(s) for s in sorted(list(open_sets), key=len)]}")

# A set is closed if its complement is open.
# In the discrete topology, every set is also closed.
closed_sets = {X - open_set for open_set in open_sets}
print(f"Closed sets (complements of open sets): {[set(s) for s in sorted(list(closed_sets), key=len)]}")

# Proper closed subsets are closed sets that are not equal to the whole space X.
proper_closed_sets = {s for s in closed_sets if s != X}
print(f"Proper closed subsets: {[set(s) for s in proper_closed_sets]}")

# Check if X can be written as a union of two proper closed subsets.
found_decomposition = False
# We search for two sets Z1 and Z2 in our list of proper_closed_sets
for Z1 in proper_closed_sets:
    for Z2 in proper_closed_sets:
        # Check if their union equals the original space X
        if Z1.union(Z2) == X:
            print("\nFound a decomposition of X into proper closed subsets:")
            
            # Format the sets for printing the equation
            x_items = sorted(list(X))
            x_str = "{" + ", ".join(map(str, x_items)) + "}"
            
            z1_items = sorted(list(Z1))
            z1_str = "{" + ", ".join(map(str, z1_items)) + "}"
            
            z2_items = sorted(list(Z2))
            z2_str = "{" + ", ".join(map(str, z2_items)) + "}"

            # Print the final equation with each number
            print(f"The equation is: {x_str} = {z1_str} U {z2_str}")
            found_decomposition = True
            break
    if found_decomposition:
        break

print("\nSince the 2-point space can be expressed as a union of its proper closed subsets, it is not irreducible.")
print("\n--- Conclusion ---")
print("We showed that 0-point and 1-point spaces are always irreducible.")
print("We constructed a 2-point space that is not irreducible.")
print("Therefore, the smallest non-negative integer n is 2.")
