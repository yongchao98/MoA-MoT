def print_filled_nilpotent_group_classification():
    """
    This function prints the classification of finite filled nilpotent groups.
    
    A finite group G is called a "filled group" if the union of all its maximal
    by inclusion product-free sets is equal to the set of its non-identity elements.
    A finite group is "nilpotent" if it is the direct product of its Sylow p-subgroups.
    """
    
    print("The classification of finite filled nilpotent groups is as follows:")
    print("-" * 70)
    
    print("A finite nilpotent group G is filled if and only if each of its Sylow p-subgroups is a filled p-group.")
    print("\nThis reduces the problem to classifying the finite filled p-groups, which are classified as follows:\n")
    
    print("Case 1: p is an odd prime")
    print("A finite p-group P for an odd prime p is filled if and only if it has an exponent of p.")
    print("(This means g^p = e for every element g in P).")
    print("Examples include the cyclic group Z_p and the elementary abelian group (Z_p)^k.\n")
    
    print("Case 2: p = 2")
    print("A finite 2-group P is filled if and only if it belongs to one of the following three classes:")
    
    print("\n  (a) P has an exponent of 2.")
    print("      (i.e., P is an elementary abelian 2-group, isomorphic to (Z_2)^k for some k >= 1).")
    print("      Example: Z_2 x Z_2 (the Klein four-group).\n")
    
    print("  (b) P is a cyclic 2-group.")
    print("      (i.e., P is isomorphic to Z_{2^k} for some k >= 1).")
    print("      Example: Z_4 or Z_8.\n")
    
    print("  (c) P is a generalized quaternion group.")
    print("      (i.e., P is isomorphic to Q_{2^n} for some n >= 3).")
    print("      Example: Q_8, the quaternion group of order 8.\n")
    
    print("-" * 70)
    print("Summary:")
    print("A finite nilpotent group G = P_1 x P_2 x ... x P_k is filled if and only if each")
    print("Sylow subgroup P_i (for prime p_i) satisfies the conditions listed above.")
    print("For example, Z_15 = Z_3 x Z_5 is filled because Z_3 and Z_5 are p-groups of exponent p.")
    print("Z_20 = Z_4 x Z_5 is filled because Z_4 is a cyclic 2-group and Z_5 is a 5-group of exponent 5.")

# Execute the function to print the classification.
if __name__ == "__main__":
    print_filled_nilpotent_group_classification()
