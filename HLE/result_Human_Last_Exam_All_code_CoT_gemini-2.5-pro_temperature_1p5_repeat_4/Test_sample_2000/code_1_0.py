def demonstrate_unbounded_ghtw(target_width):
    """
    Demonstrates that the generalized hypertreewidth for a hypergraph
    with 3 edges and unbounded rank is itself unbounded.

    Args:
        target_width (int): The target generalized hypertreewidth to demonstrate.
    """
    if not isinstance(target_width, int) or target_width < 0:
        print("Error: Please provide a non-negative integer for the target width.")
        return

    # To construct a hypergraph with a generalized hypertreewidth of 'target_width',
    # we need its largest hyperedge to have a size of 'target_width + 1'.
    e1_size = target_width + 1

    # Define the three hyperedges. Let e1 be the large hyperedge.
    # The vertices can be represented by integers {0, 1, ..., target_width}.
    e1 = set(range(e1_size))
    # The other two hyperedges can be empty or subsets of e1.
    e2 = set()
    e3 = set()

    hypergraph = [e1, e2, e3]

    print(f"Let's construct a hypergraph H with 3 hyperedges for a target ghtw of {target_width}:")
    # To keep the output clean, we'll only show a few elements of e1 if it's large.
    e1_str = str(e1)
    if len(e1_str) > 60:
        e1_str = f"{{0, 1, 2, ..., {e1_size-1}}}"
    print(f"e1 = {e1_str}")
    print(f"e2 = {e2}")
    print(f"e3 = {e3}")
    print("-" * 30)

    print("Reasoning for the hypertreewidth of this hypergraph:")
    print("1. Any valid generalized tree decomposition (GHTD) must have at least one bag, say B_t,")
    print(f"   such that e1 is a subset of B_t. This is required by the definition of a GHTD.")
    print(f"2. Therefore, the size of this bag, |B_t|, must be at least the size of e1.")
    print(f"   |B_t| >= |e1| = {e1_size}.")
    print("3. The width of any GHTD is defined as (max |Bag|) - 1.")
    print(f"   This means the width must be at least {e1_size} - 1 = {target_width}.")
    print(f"   So, the generalized hypertreewidth, ghtw(H), is >= {target_width}.")
    print("4. We can achieve this lower bound by using a simple GHTD with a single bag B_0 = e1.")
    print("   This decomposition is valid and has a width of |e1| - 1.")
    print("-" * 30)
    
    print("Conclusion:")
    print(f"The generalized hypertreewidth of our constructed hypergraph is exactly {target_width}.")
    print("Since we can choose any non-negative integer for 'target_width', the maximum")
    print("generalized hypertreewidth for a hypergraph with 3 edges is unbounded.")
    print("-" * 30)
    
    # As requested, output the final equation for the example calculation.
    print("The final equation for the ghtw of the constructed example:")
    print(f"{target_width} = {e1_size} - 1")


# You can change this value to see the construction for any target width.
example_target_width = 2
demonstrate_unbounded_ghtw(example_target_width)
