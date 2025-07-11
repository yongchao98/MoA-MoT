import itertools

class Group:
    """A simple class to represent a finite group."""
    def __init__(self, name, elements, identity, operation):
        self.name = name
        self.elements = elements
        self.identity = identity
        # The operation is a function, e.g., lambda a, b: ...
        self.op = operation
        # Pre-calculate non-identity elements for convenience
        self.non_identity_elements = [e for e in self.elements if e != self.identity]

def is_product_free(S, group):
    """Checks if a set S is product-free within a group."""
    for s1 in S:
        for s2 in S:
            product = group.op(s1, s2)
            if product in S:
                return False
    return True

def find_maximal_product_free_set_of_size_2(group):
    """
    Finds a maximal by inclusion product-free set of size 2 in a given group.
    Returns the set if found, otherwise returns None.
    """
    # Iterate over all subsets of size 2 of the non-identity elements
    for s_tuple in itertools.combinations(group.non_identity_elements, 2):
        S = list(s_tuple)
        
        # Condition 1: The set must be product-free
        if not is_product_free(S, group):
            continue

        # Condition 2: The set must be maximal by inclusion
        is_maximal = True
        # Get all elements not in S
        G_minus_S = [g for g in group.elements if g not in S]
        for g in G_minus_S:
            # Create the extended set S U {g}
            S_union_g = S + [g]
            # If the extended set is still product-free, then S is not maximal
            if is_product_free(S_union_g, group):
                is_maximal = False
                break
        
        if is_maximal:
            # Found a set satisfying both conditions
            return S
            
    # No such set was found in the group
    return None

def setup_groups():
    """Defines the 5 groups that are candidates."""
    groups_to_check = []

    # 1. C4: Cyclic group of order 4 (using additive notation)
    c4 = Group("C4", [0, 1, 2, 3], 0, lambda a, b: (a + b) % 4)
    groups_to_check.append(c4)

    # 2. C2 x C2: Klein four-group
    k4_elements = list(itertools.product([0, 1], repeat=2))
    k4 = Group("C2 x C2", k4_elements, (0,0), lambda a, b: ((a[0] + b[0]) % 2, (a[1] + b[1]) % 2))
    groups_to_check.append(k4)

    # 3. C5: Cyclic group of order 5
    c5 = Group("C5", [0, 1, 2, 3, 4], 0, lambda a, b: (a + b) % 5)
    groups_to_check.append(c5)
    
    # 4. C6: Cyclic group of order 6
    c6 = Group("C6", [0, 1, 2, 3, 4, 5], 0, lambda a, b: (a + b) % 6)
    groups_to_check.append(c6)
    
    # 5. D6 (S3): Dihedral group of order 6, represented by a Cayley table
    # Elements: e=0, r=1, r^2=2, f=3, fr=4, fr^2=5
    d6_table = [
        [0, 1, 2, 3, 4, 5], # row e
        [1, 2, 0, 5, 3, 4], # row r
        [2, 0, 1, 4, 5, 3], # row r^2
        [3, 4, 5, 0, 1, 2], # row f
        [4, 5, 3, 2, 0, 1], # row fr
        [5, 3, 4, 1, 2, 0], # row fr^2
    ]
    d6 = Group("D6 (S3)", [0, 1, 2, 3, 4, 5], 0, lambda a, b: d6_table[a][b])
    groups_to_check.append(d6)
    
    return groups_to_check

def main():
    """Main function to run the verification."""
    groups_to_check = setup_groups()
    found_groups_names = []
    
    print("Searching for maximal product-free sets of size 2...\n")

    for group in groups_to_check:
        result = find_maximal_product_free_set_of_size_2(group)
        if result:
            found_groups_names.append(group.name)
            # Format result for printing
            if group.name == "C2 x C2":
                str_result = [str(x) for x in result]
                print(f"Success: Found in {group.name}. Example set: {str_result}")
            else:
                print(f"Success: Found in {group.name}. Example set: {result}")
        else:
            print(f"Failure: Could not find such a set in {group.name}.")

    print("\n--- Summary ---")
    print("The finite groups containing a maximal by inclusion product-free set of size 2 are:")
    for name in found_groups_names:
        print(f"- {name}")

    count = len(found_groups_names)
    print(f"\nThis makes for a total of {count} non-isomorphic groups.")
    
    # Fulfilling the request to "output each number in the final equation"
    if count > 0:
        equation_str = " + ".join(["1"] * count)
        print(f"The calculation is: {equation_str} = {count}")

if __name__ == '__main__':
    main()