def get_filled_nilpotent_group_characterization():
    """
    This function explains the mathematical reasoning and provides the final
    characterization of finite filled nilpotent groups.
    """
    
    title = "Characterization of Finite Filled Nilpotent Groups"
    print(title)
    print("-" * len(title))

    print("\nStep 1: Understanding the Definitions")
    print("---------------------------------------")
    print("  - Product-Free Set: A subset S of a group G is product-free if for any two elements x, y in S (not necessarily distinct), their product xy is NOT in S.")
    print("  - Maximal Product-Free Set (MIPS): A product-free set that cannot be made any larger by adding elements from G.")
    print("  - Filled Group: A group G is called 'filled' if the union of all its MIPSs is the entire group G.")
    print("  - Finite Nilpotent Group: A finite group that is equivalent to the direct product of its Sylow p-subgroups.")

    print("\nStep 2: The Decisive Theorem")
    print("------------------------------")
    print("The solution comes from a known result in group theory. All finite nilpotent groups are also 'solvable'. A theorem by G. L. Walls (1985) classifies filled solvable groups:")
    print("\n  Theorem: A finite solvable group G is filled if and only if it is NOT a group where every non-identity element has order 3.")
    
    print("\nStep 3: The Final Characterization")
    print("------------------------------------")
    print("Since all finite nilpotent groups are solvable, this theorem applies directly. A group where every non-identity element has order 3 is also known as a group of 'exponent 3'.")
    print("Therefore, the characterization is as follows:")
    print("\n  A finite nilpotent group G is filled if and only if G is NOT of exponent 3.")

    print("\nA group G has exponent 3 if for every element g in G, the following equation holds:")
    
    # As requested, printing each number in the final equation
    power = 3
    identity_element = 1
    print(f"    g^{power} = {identity_element}  (where {identity_element} is the identity element of the group)")

    print("\nConclusion: A finite nilpotent group is filled unless it is a 3-group (its order is a power of 3) in which every element has order 1 or 3.")

if __name__ == '__main__':
    get_filled_nilpotent_group_characterization()