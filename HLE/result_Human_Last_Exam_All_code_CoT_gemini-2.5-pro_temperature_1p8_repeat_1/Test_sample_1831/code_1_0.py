import math

def main():
    """
    Calculates the number of equivalence classes of endomorphisms on a set of size 4.
    This is equivalent to counting the number of non-isomorphic functional graphs on 4 vertices.
    The method is to sum the counts over all integer partitions of 4.
    """
    print("This script calculates the number of non-isomorphic functional graphs on 4 vertices.")
    print("The approach is to analyze the structure of these graphs based on integer partitions of 4.\n")

    # Helper values: Number of non-isomorphic rooted trees on n vertices.
    # A connected functional graph with a 1-cycle is just a rooted tree.
    # t_n = number of rooted trees on n vertices.
    num_rooted_trees = {
        1: 1,  # A single node.
        2: 1,  # A root with one child.
        3: 2,  # A path (length 2) and a root with two children.
        4: 4   # See OEIS A000081.
    }

    # Memoization for the number of connected functional graphs on n vertices.
    memo_graphs = {}

    def count_connected_graphs(n):
        """
        Counts the number of non-isomorphic CONNECTED functional graphs on n vertices.
        A connected functional graph has a single cycle.
        """
        if n in memo_graphs:
            return memo_graphs[n]

        if n == 0:
            return 1 # Base case for calculations
        if n == 1:
            # Must be a 1-cycle (fixed point).
            return 1

        total = 0
        # Sum over possible cycle lengths k from 1 to n.
        # Case k=1: The graph is a rooted tree.
        total += num_rooted_trees[n]
        if n > 1:
            # Case k=2: 2-cycle with n-2 vertices in trees attached.
            # Partition n-2 vertices between 2 cycle nodes.
            # - (n-2) + 0: forms a rooted tree of size n-1. Num ways = t_{n-1}.
            # - i + (n-2-i): smaller partitions.
            if n == 2: total += 1 # 2-cycle
            if n == 3: total += num_rooted_trees[2] # t_2=1
            if n == 4: total += num_rooted_trees[3] + 1 # t_3 + (1+1 ways) = 2+1=3
        if n > 2:
            # Case k=3: 3-cycle with n-3 vertices attached.
            if n == 3: total += 1 # 3-cycle
            if n == 4: total += num_rooted_trees[2] # t_2=1 way to attach 1 vertex
        if n > 3:
            # Case k=4: 4-cycle
            if n == 4: total += 1

        memo_graphs[n] = total
        return total

    # --- Calculations for each partition of 4 ---

    total_classes = 0

    # Partition [4]: One component of size 4.
    count_part4 = count_connected_graphs(4)
    print(f"Partition [4]: Number of connected graphs on 4 vertices = {count_part4}")
    total_classes += count_part4

    # Partition [3, 1]: One component of size 3 and one of size 1.
    count_part3_1 = count_connected_graphs(3) * count_connected_graphs(1)
    print(f"Partition [3, 1]: (Graphs on 3 vertices) * (Graphs on 1 vertex) = {count_connected_graphs(3)} * {count_connected_graphs(1)} = {count_part3_1}")
    total_classes += count_part3_1

    # Partition [2, 2]: Two components of size 2.
    num_types_2 = count_connected_graphs(2)
    # This is combinations with repetition: choose 2 items from `num_types_2` types.
    # C(n+k-1, k) where n=num_types_2, k=2.
    count_part2_2 = math.comb(num_types_2 + 2 - 1, 2)
    print(f"Partition [2, 2]: Ways to choose 2 components from {num_types_2} types of size-2 graphs = {count_part2_2}")
    total_classes += count_part2_2

    # Partition [2, 1, 1]: One component of size 2, two of size 1.
    # The two size-1 components are identical (fixed points).
    # We just choose the type for the size-2 component.
    count_part2_1_1 = count_connected_graphs(2)
    print(f"Partition [2, 1, 1]: Number of choices for the size-2 component = {count_part2_1_1}")
    total_classes += count_part2_1_1

    # Partition [1, 1, 1, 1]: Four components of size 1.
    # Only one way: four fixed points.
    count_part1_1_1_1 = 1
    print(f"Partition [1, 1, 1, 1]: Only {count_part1_1_1_1} way (four fixed points)")
    total_classes += count_part1_1_1_1

    # --- Final Result ---
    print("\n----------------------------------------------------")
    print("The total number of equivalence classes is the sum of the counts for each partition.")
    print(f"Final calculation: {count_part4} + {count_part3_1} + {count_part2_2} + {count_part2_1_1} + {count_part1_1_1_1} = {total_classes}")
    print("----------------------------------------------------")


if __name__ == "__main__":
    main()
<<<19>>>