import sys

def print_obstruction_groups(n_symbol='n', k_symbol='k'):
    """
    This function prints the list of homotopy and homology groups that classify
    the obstructions for two paths of bundle automorphisms to be homotopic.

    Args:
      n_symbol (str): The symbol for the dimension of the base manifold's suspension.
      k_symbol (str): The symbol for the rank of the vector bundle.
    """
    # Validate that k is a positive integer for the problem to be well-defined.
    # The actual value is not used, only its symbol.
    try:
        rank = 2 * int(k_symbol)
    except ValueError:
        rank_symbol = f"2*{k_symbol}"

    n_minus_1 = f"{n_symbol}-1"
    
    print("The obstructions for the two paths to be homotopic are classified by a list of algebraic topology groups.")
    print("These groups arise from two sources:")
    print("\n1. A 'fiberwise' obstruction, which is constant across the base space. It is classified by the group:")
    print(f"   \u03c0\u2081(SO({rank_symbol}))")
    print("\n2. A 'twisting' obstruction, which captures how the fiberwise loops are twisted over the base space \u03A3X. This obstruction is described by groups built from the homology of X and the homotopy of SO(2k):")
    
    group_list_title = f"   H^{n_minus_1}(X; \u03c0\u209a(SO({rank_symbol}))) for q \u2265 1"
    print(group_list_title)
    
    print("\nFor example, the first few groups in this family are:")
    for q in range(1, 5):
      print(f"   H^({n_minus_1})(X; \u03c0\u208{q}(SO({rank_symbol})))")
    print("   ...")
    print("\nNote: For n > 2, since X is a homology sphere, H^{n-1}(X; G) is isomorphic to G for any coefficient group G. In this case, the second list of groups simplifies to \u03c0\u209a(SO(2k)) for q \u2265 1.")


if __name__ == '__main__':
    # The problem is given in terms of symbolic variables n and k.
    # We will call the function with these symbols.
    print_obstruction_groups('n', 'k')
