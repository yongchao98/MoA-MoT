def describe_filled_groups():
    """
    This program describes the nonabelian filled groups of order 2*q^m,
    where q is an odd prime and m is a natural number, based on known
    classification theorems in group theory.
    """

    # We use symbolic representations for the odd prime 'q' and natural number 'm'.
    q_str = 'q'
    m_str = 'm'

    print("The nonabelian filled groups G of order |G| = 2 * q^m, for an odd prime q and natural number m,")
    print("are the dihedral groups of that order.")
    print("-" * 70)

    # Define the parameters for the dihedral group D_2n
    n_str = f"{q_str}^{m_str}"
    order_str = f"2 * ({n_str})"

    print(f"These groups are the dihedral groups, denoted D_({order_str}).")
    print(f"Here, n = {n_str}, and the subscript {order_str} represents the order of the group.")
    print("\nThis group can be defined by a presentation, which consists of generators and relations.")

    # The presentation
    print("\nGroup Presentation:")
    print(f"Let n = {n_str}. The presentation is <r, s | r^n = e, s^2 = e, srs = r^(-1)>")
    print("where 'e' is the identity element, 'r' represents a rotation, and 's' represents a reflection.")

    print("\nThe defining relations (the equations) are:")

    # Relation 1: r^n = e
    r_exponent = n_str
    print(f"\n1. The rotation element 'r' has order n = {n_str}:")
    print(f"   r^({r_exponent}) = e")

    # Relation 2: s^2 = e. This prints the number '2'.
    s_exponent = 2
    print(f"\n2. The reflection element 's' has order {s_exponent}:")
    print(f"   s^{s_exponent} = e")

    # Relation 3: srs = r^-1. This shows the number '-1' as an exponent.
    print("\n3. The conjugation of 'r' by 's' gives the inverse of 'r':")
    print(f"   s * r * s = r^(-1)")

    print("-" * 70)
    print("Summary:")
    print(f"For any odd prime q and any natural number m, the only nonabelian filled group")
    print(f"of order 2*q^m is the dihedral group D_(2*q^m).")
    print("\nFor a concrete example, if q=3 and m=1, the group is D_6 (order 6),")
    print("which is isomorphic to the symmetric group S_3.")


if __name__ == '__main__':
    describe_filled_groups()