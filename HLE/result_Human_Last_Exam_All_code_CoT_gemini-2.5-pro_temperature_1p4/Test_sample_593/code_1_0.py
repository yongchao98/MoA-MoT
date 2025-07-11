def calculate_treewidth_bound(t_H, t_G, k):
    """
    Calculates the tight upper bound on the treewidth of the union of two graphs.

    Let H be a graph of treewidth t_H, and G be a graph of treewidth t_G.
    Assume they share k vertices, S = V(H) intersect V(G).
    The graph F is the union of H and G.
    A tight upper bound on the treewidth of F, t_F, is max(t_H, t_G) + k.

    Args:
        t_H (int): The treewidth of graph H.
        t_G (int): The treewidth of graph G.
        k (int): The number of vertices in the intersection of H and G.
    """
    if not all(isinstance(i, int) and i >= 0 for i in [t_H, t_G, k]):
        print("Treewidth and k must be non-negative integers.")
        return

    # The formula for the tight upper bound
    max_t = max(t_H, t_G)
    bound = max_t + k

    print(f"Given values: t_H = {t_H}, t_G = {t_G}, k = {k}")
    print("The tight upper bound for the treewidth of the combined graph F is given by the formula:")
    print("tw(F) <= max(t_H, t_G) + k")
    print("\nCalculation:")
    
    # Printing the equation with substituted values
    # f"max({t_H}, {t_G}) + {k} = {max_t} + {k} = {bound}" also works, but the below format is clearer.
    final_equation_step1 = f"max({t_H}, {t_G}) + {k}"
    final_equation_step2 = f"{max_t} + {k}"
    final_equation_step3 = f"{bound}"
    
    print(f"{final_equation_step1} = {final_equation_step2}")
    print(f"{final_equation_step2} = {final_equation_step3}")
    print(f"\nThe result is: {bound}")


# Example usage with some values.
# You can change these values to see the result for different graphs.
example_t_H = 7
example_t_G = 4
example_k = 3

calculate_treewidth_bound(example_t_H, example_t_G, example_k)
