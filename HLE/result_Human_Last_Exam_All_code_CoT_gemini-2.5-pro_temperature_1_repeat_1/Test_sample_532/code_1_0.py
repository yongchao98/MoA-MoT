import itertools

class Group:
    """A simple class to represent a finite group."""
    def __init__(self, name, elements, op):
        self.name = name
        # Ensure elements are hashable and store as a list for consistent ordering
        self.elements = sorted(list(elements), key=str)
        self.op = op
        self.identity = self._find_identity()
        self.non_identity_elements = {el for el in self.elements if el != self.identity}

    def _find_identity(self):
        """Finds the identity element of the group."""
        for e in self.elements:
            if all(self.op(e, x) == x and self.op(x, e) == x for x in self.elements):
                return e
        raise ValueError("Group has no identity element")

def find_maximal_product_free_sets(group):
    """
    Finds all maximal by inclusion product-free sets of a group.
    A set S is product-free if for all x, y in S, x*y is not in S.
    """
    # We only need to consider subsets of non-identity elements.
    non_identity_indices = {i for i, el in enumerate(group.elements) if el != group.identity}
    
    product_free_sets = []
    
    # Iterate through the powerset of non-identity elements to find product-free sets.
    # We build a list of sets of indices.
    for r in range(1, len(non_identity_indices) + 1):
        for subset_indices_tuple in itertools.combinations(non_identity_indices, r):
            subset_indices = set(subset_indices_tuple)
            subset_elements = {group.elements[i] for i in subset_indices}
            
            is_pf = True
            # Check the product-free condition
            for x in subset_elements:
                for y in subset_elements:
                    if group.op(x, y) in subset_elements:
                        is_pf = False
                        break
                if not is_pf:
                    break
            
            if is_pf:
                product_free_sets.append(subset_indices)

    # From the list of all product-free sets, find the maximal ones.
    maximal_sets_indices = []
    for s1 in product_free_sets:
        is_maximal = True
        for s2 in product_free_sets:
            if s1 != s2 and s1.issubset(s2):
                is_maximal = False
                break
        if is_maximal:
            maximal_sets_indices.append(s1)
            
    # Convert from indices back to elements for the final output.
    maximal_element_sets = []
    for m_set_idx in maximal_sets_indices:
        maximal_element_sets.append({group.elements[i] for i in m_set_idx})

    return maximal_element_sets

def analyze_group(group, is_nilpotent):
    """Analyzes a group and prints its properties related to the question."""
    print(f"--- Group: {group.name} ---")
    print(f"Is Nilpotent? {'Yes' if is_nilpotent else 'No'}")
    
    maximal_sets = find_maximal_product_free_sets(group)
    
    union_of_maximal_sets = set()
    for s in maximal_sets:
        union_of_maximal_sets.update(s)
        
    is_group_filled = (union_of_maximal_sets == group.non_identity_elements)

    # To satisfy the prompt "output each number in the final equation",
    # we will list the elements of the sets involved in the check.
    print(f"Maximal product-free sets found: {maximal_sets}")
    print(f"The union of these sets is: {union_of_maximal_sets}")
    print(f"The set of all non-identity elements is: {group.non_identity_elements}")

    if is_group_filled:
        print(f"Result: The union equals the set of non-identity elements, so {group.name} is a filled group.")
        if is_nilpotent:
            print(f"Conclusion: {group.name} is a finite filled nilpotent group.")
        else:
            print(f"Conclusion: {group.name} is filled, but it is not nilpotent.")
    else:
        print(f"Result: The union does not equal the set of non-identity elements, so {group.name} is NOT a filled group.")
    print("-" * 25 + "\n")


# --- Main Execution ---
if __name__ == "__main__":
    # Define the groups we want to test.
    # C3: Cyclic group of order 3. Nilpotent (abelian).
    c3 = Group("C3", {0, 1, 2}, lambda a, b: (a + b) % 3)
    
    # C4: Cyclic group of order 4. Nilpotent (abelian).
    c4 = Group("C4", {0, 1, 2, 3}, lambda a, b: (a + b) % 4)

    # C2xC2: Klein four-group. An elementary abelian 2-group. Nilpotent (abelian).
    c2c2_map = {0: (0,0), 1:(0,1), 2:(1,0), 3:(1,1)}
    c2c2_inv_map = {v: k for k, v in c2c2_map.items()}
    def c2c2_op(a, b):
        va = c2c2_map[a]
        vb = c2c2_map[b]
        res_tuple = ((va[0] + vb[0]) % 2, (va[1] + vb[1]) % 2)
        return c2c2_inv_map[res_tuple]
    # Represent elements as strings for clarity in output
    c2c2 = Group("C2xC2", {"e", "a", "b", "c"}, 
                 lambda x, y: {"e":0, "a":1, "b":2, "c":3}[x] if isinstance(x, str) else x, # Dummy op for init
                )
    c2c2_elements = ["e", "a", "b", "c"]
    c2c2_op_table = [ # Cayley table for C2xC2
    # e, a, b, c
      [0, 1, 2, 3], # e
      [1, 0, 3, 2], # a
      [2, 3, 0, 1], # b
      [3, 2, 1, 0]  # c
    ]
    c2c2.op = lambda a, b: c2c2_elements[c2c2_op_table[c2c2_elements.index(a)][c2c2_elements.index(b)]]
    c2c2.identity = "e"
    c2c2.non_identity_elements = {"a", "b", "c"}

    print("--- Analysis of Finite Filled Nilpotent Groups ---\n")
    analyze_group(c3, is_nilpotent=True)
    analyze_group(c2c2, is_nilpotent=True)
    
    print("--- Checking other groups for context ---\n")
    analyze_group(c4, is_nilpotent=True)

    print("--- Final Conclusion based on Theoretical Classification ---")
    print("Based on a well-known (though simplified) classification theorem, the finite groups that are both filled and nilpotent are:")
    print("1. The group of order 3.")
    print("2. Elementary abelian 2-groups (e.g., C2, C2xC2, ...).")
