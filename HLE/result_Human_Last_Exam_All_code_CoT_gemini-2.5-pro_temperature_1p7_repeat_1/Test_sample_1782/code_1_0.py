import random

class Partition:
    """Represents a partition of a set of integers."""
    def __init__(self, sets):
        # A partition is a list of sets.
        self.parts = [frozenset(p) for p in sets]

    def __str__(self):
        # String representation for printing.
        # To avoid clutter, we'll only show the number of sets in the partition.
        return f"A partition with {len(self.parts)} sets."

    def get_all_elements(self):
        """Returns the union of all sets in the partition."""
        return set().union(*self.parts)

def refine_partition(partition):
    """
    Refines a given partition by splitting each of its sets into two non-empty subsets.
    """
    new_parts = []
    for part in partition.parts:
        if len(part) < 2:
            # Cannot split a set with one element
            new_parts.append(part)
            continue
        
        # Convert to list to be able to draw elements
        part_list = list(part)
        random.shuffle(part_list)
        
        # Split point is random, but ensuring it's not at the ends
        # to guarantee a non-trivial split.
        split_point = random.randint(1, len(part_list) - 1)
        
        part1 = set(part_list[:split_point])
        part2 = set(part_list[split_point:])
        
        new_parts.append(part1)
        new_parts.append(part2)
        
    return Partition(new_parts)

def main():
    """
    Main function to simulate the refinement process.
    """
    N = 100  # Our stand-in for "omega_1"
    base_set = set(range(N))
    
    # L_0: Start with a simple partition of the base set into two halves.
    print("This script simulates building a tree of refining partitions on a FINITE set.")
    print("This is an analogy, not a representation, of the situation in P(omega_1)/<omega_1>.")
    print("-" * 20)

    initial_partition = Partition([set(range(N//2)), set(range(N//2, N))])
    
    # The "tree" is a sequence of partitions (levels).
    tree_levels = [initial_partition]
    
    # We create a sequence of 10 refinements.
    num_levels = 10
    for i in range(num_levels):
        new_level = refine_partition(tree_levels[-1])
        tree_levels.append(new_level)

    # Print the description of the generated levels
    for i, level in enumerate(tree_levels):
        print(f"Level L_{i}: {level}")

    print("-" * 20)
    print("Why this finite analogy has a common refinement:")
    print("In this finite model, the refinement process must terminate when all sets in the partition")
    print("are singletons (sets with one element). The collection of all these singleton sets,")
    print("e.g., {{0}, {1}, {2}, ..., {N-1}}, forms a partition that refines every level")
    print("of the tree we constructed. This is called a 'common refinement'.")
    print("\nThe crucial result from set theory is that for the uncountable set omega_1, it's")
    print("possible to construct a similar tower of `omega_1` refining partitions in such a clever")
    print("way that NO such common refinement can exist. This is a profound property of uncountable sets.")

if __name__ == "__main__":
    main()
