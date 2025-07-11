def explain_order_type():
    """
    Explains the derivation of the ordinal equation for the order type
    of lexically ordered finite strings on {a,b,c,d}.
    """
    k = 4
    alphabet = "{a, b, c, d}"

    print("Let Omega be the order type of the set S of finite strings on {}.".format(alphabet))
    print("The set S can be partitioned into {} subsets based on the first character.".format(k))
    print("S = (Strings starting with 'a') U (Strings starting with 'b') U ...")
    print("Lexicographical order means these are ordered blocks.")
    print("So, Omega = Type(S_a) + Type(S_b) + Type(S_c) + Type(S_d).")
    print("-" * 20)
    print("Consider the structure of one subset, S_a:")
    print("S_a = {'a'} U {'a' + s for s in S}")
    print("The first element 'a' has order type 1.")
    print("The set {'a' + s | s in S} is order-isomorphic to S, so its type is Omega.")
    print("Thus, Type(S_a) = 1 + Omega.")
    print("-" * 20)
    print("By symmetry, all {} partitions have type (1 + Omega).".format(k))
    print("Substituting this back gives the final ordinal equation:")
    print("Omega = (1 + Omega) + (1 + Omega) + (1 + Omega) + (1 + Omega)")
    print("Which can be written as:")
    
    one = 1
    four = 4
    # The final equation printed for the user, as requested.
    print(f"Omega = ({one} + Omega) * {four}")
    
    print("\nThe smallest transfinite ordinal solving this equation is omega^omega.")

explain_order_type()