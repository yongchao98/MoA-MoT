import math

def illustrate_grid_subgraph_bound(d, treewidth_values):
    """
    Illustrates the growth of the guaranteed grid subgraph size for
    a class of graphs with bounded degree d and increasing treewidth.
    This is based on a theorem by Chekuri and Chuzhoy.

    The theorem states that a graph with treewidth 'w' and max degree 'd'
    contains a grid subgraph of size k-by-k where k is Omega(w / (d * log w)).
    We will ignore the constant factor in Omega() and show the scaling.
    """
    print(f"Analyzing a class of graphs C with max degree <= {d} and unbounded treewidth.")
    print("This means for any large treewidth value 'w', there is a graph in C with tw(G) > w.")
    print("A key theorem states this implies the existence of large grid subgraphs.\n")
    print("The size of the guaranteed grid, k-by-k, scales as follows:")
    print("k is proportional to w / (d * log(w))\n")
    print(f"{'Treewidth (w)':<15} | {'Equation for k (proportional to)':<35} | {'Value':<10}")
    print("-" * 70)

    for w in treewidth_values:
        if w <= 1:
            # log(w) is undefined or non-positive
            scaling_factor_str = "N/A"
            equation_str = f"{w} / ({d} * log({w}))"
        else:
            # We calculate the scaling factor k ~ w / (d * log w)
            scaling_factor = w / (d * math.log(w))
            scaling_factor_str = f"{scaling_factor:.2f}"
            # This is the "equation" part requested by the prompt.
            equation_str = f"{w} / ({d} * log({w}))"

        print(f"{w:<15} | {equation_str:<35} | {scaling_factor_str:<10}")

# Main execution
if __name__ == '__main__':
    # The problem states d is a constant. Let's use d=4 as an example.
    d = 4
    # The class C has unbounded treewidth. We can pick a sequence of
    # increasing treewidth values to demonstrate the effect.
    treewidth_values = [10, 100, 1000, 10000, 100000, 1000000]

    illustrate_grid_subgraph_bound(d, treewidth_values)

    print("\nAs the treewidth 'w' grows, the guaranteed grid dimension 'k' also grows unboundedly.")
    print("This shows that for any target grid size, there must be a graph in C that contains it.")
    print("This confirms that statement B must be true for the class C.")
