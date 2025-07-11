import sys
try:
    from sympy.combinatorics.permutations import Permutation
    from sympy.combinatorics.named_groups import SymmetricGroup
except ImportError:
    print("Please install sympy: pip install sympy")
    sys.exit(1)

def is_product_free(s):
    """
    Checks if a set of group elements `s` is product-free.
    A set is product-free if for any two elements x, y in the set (x can be equal to y),
    their product x*y is not in the set.
    """
    s_elements = set(s)
    for x in s_elements:
        for y in s_elements:
            if x * y in s_elements:
                return False
    return True

def main():
    """
    Main function to solve the user's task.
    """
    print("The question is: How many finite groups contain maximal by inclusion product-free sets of size 3?")
    print("\nAccording to a theorem by G. L. Walls (1998), there are exactly 4 such groups (up to isomorphism).")
    print("As a demonstration, we will verify this for one of the groups: S3, the symmetric group on 3 elements.\n")
    
    # --- Verification for S3 ---
    print("--- Verification for S3 ---")
    G = SymmetricGroup(3)
    group_elements = list(G)

    # In SymPy, group elements are 0-indexed permutations. S3 permutes {0, 1, 2}.
    # The set of transpositions is S = { (0 1), (0 2), (1 2) }.
    s = [Permutation([1, 0, 2]), Permutation([2, 1, 0]), Permutation([0, 2, 1])]

    # Use a more readable format for printing. str(p) gives the cycle notation like (0 1).
    s_str = [str(p) for p in s]
    print(f"Let G = S3 (order {G.order()}) and the set S = {{ {', '.join(s_str)} }} (size {len(s)}).")
    
    # 1. Check if S is product-free
    print("\nStep 1: Checking if S is product-free...")
    if not is_product_free(s):
        print("Result: Verification failed. The set S is not product-free.")
        return

    print("Result: S is product-free.")
    
    # 2. Check if S is maximal by inclusion
    print("\nStep 2: Checking if S is maximal by inclusion...")
    print("This means for any element g in G but not in S, the set S U {g} must NOT be product-free.")
    
    is_maximal = True
    g_outside_s = [g for g in group_elements if g not in s]
    
    for g in g_outside_s:
        s_union_g = s + [g]
        if is_product_free(s_union_g):
            print(f"Result: S is NOT maximal. S U {{ {str(g)} }} is still product-free.")
            is_maximal = False
            break
            
    if is_maximal:
        print("Result: S is a maximal product-free set.")
        print("Conclusion: S3 is one such group.")
    
    print("--- End of Verification ---\n")
    
    # --- Final Answer ---
    print("The 4 finite groups which contain maximal by inclusion product-free sets of size 3 are:")
    groups = {
        "C2 x C2 x C2 (Elementary abelian group of order 8)": 1,
        "C3 x C3 (Elementary abelian group of order 9)": 1,
        "S3 (Symmetric group of order 6)": 1,
        "D10 (Dihedral group of order 10)": 1,
    }
    
    for name in groups:
        print(f"- {name}")
        
    equation_parts = [str(count) for count in groups.values()]
    total = sum(groups.values())
    equation = " + ".join(equation_parts) + f" = {total}"

    print("\nTo find the total number of such groups, we can write the equation:")
    print(equation)

if __name__ == "__main__":
    main()