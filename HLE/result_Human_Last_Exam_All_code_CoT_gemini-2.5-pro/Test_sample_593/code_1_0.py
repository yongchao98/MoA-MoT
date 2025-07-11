import sys

def solve_treewidth_union(t_H, t_G, k):
    """
    Calculates the tight upper bound on the treewidth of the union of two graphs.

    Args:
        t_H (int): The treewidth of graph H.
        t_G (int): The treewidth of graph G.
        k (int): The number of vertices in the intersection of H and G.
    """
    # The problem statement implies k >= 0. Treewidth is also >= 0,
    # except for the empty graph. We assume non-empty graphs.
    # The formula for the tight upper bound on the treewidth of F is max(t_H, t_G, k-1).
    if t_H < 0 or t_G < 0 or k < 0:
      print("Error: Treewidth and the number of shared vertices cannot be negative.", file=sys.stderr)
      return

    # For k=0, the intersection is empty. The treewidth of the disjoint union is max(t_H, t_G).
    # The term k-1 becomes -1, which is correctly handled by the max function.
    bound = max(t_H, t_G, k - 1)

    # Output the explanation and the final equation.
    print("The tight upper bound on the treewidth of F is given by the formula: t_F <= max(t_H, t_G, k - 1)")
    print("\nFor the given values:")
    print(f"t_H = {t_H}")
    print(f"t_G = {t_G}")
    print(f"k = |V(H) intersect V(G)| = {k}")
    print("\nThe calculation is:")
    # The final print statement shows the full equation as requested.
    print(f"t_F <= max({t_H}, {t_G}, {k} - 1) = {bound}")


# Example values to demonstrate the calculation.
# These values are inspired by the tightness proof where t_H and t_G are small,
# but k is large, making the bound equal to k-1.
example_t_H = 1
example_t_G = 1
example_k = 10

solve_treewidth_union(example_t_H, example_t_G, example_k)
