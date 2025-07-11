import sys

def solve_and_explain():
    """
    Finds the smallest non-negative integer n for which an n-point 
    topological space that is not irreducible exists.

    The reasoning is as follows:
    A topological space X is NOT irreducible if it can be written as a union of
    two of its proper closed subsets, say Z1 and Z2, such that X = Z1 U Z2.

    - For n=0 (X={}) and n=1 (X={p}), it can be shown that no such pair of
      proper closed subsets exists. All 0-point and 1-point spaces are irreducible.

    - For n=2, we can construct a counterexample. Let the space be X = {0, 1}
      with the discrete topology (where every subset is open, and thus also closed).
      We can choose two proper closed subsets Z1 = {0} and Z2 = {1}. Their union
      is {0} U {1} = {0, 1}, which is the entire space X.

    This proves that for n=2, a non-irreducible space exists. Therefore, 2 is the
    smallest such non-negative integer. The code below illustrates this finding.
    """
    n = 2
    
    # Define the 2-point space
    X = {0, 1}
    
    # In the discrete topology on X, {0} and {1} are proper closed subsets.
    Z1 = {0}
    Z2 = {1}
    
    # Calculate the union of these two proper closed subsets.
    union_Z1_Z2 = Z1.union(Z2)
    
    # The final equation shows that X is a union of two proper closed subsets.
    print(f"Let n = {n}. Consider the space X = {sorted(list(X))} with the discrete topology.")
    print("In this topology, every subset is closed.")
    
    z1_str = str(sorted(list(Z1))).replace('[', '{').replace(']', '}')
    z2_str = str(sorted(list(Z2))).replace('[', '{').replace(']', '}')
    union_str = str(sorted(list(union_Z1_Z2))).replace('[', '{').replace(']', '}')
    
    print(f"We can choose the proper closed subsets Z1 = {z1_str} and Z2 = {z2_str}.")
    print("\nTheir union covers the entire space X:")
    print(f"{z1_str} U {z2_str} = {union_str}")
    
    print(f"\nSince a non-irreducible space exists for n={n}, and not for n < {n},")
    print(f"the smallest such non-negative integer is {n}.")

solve_and_explain()