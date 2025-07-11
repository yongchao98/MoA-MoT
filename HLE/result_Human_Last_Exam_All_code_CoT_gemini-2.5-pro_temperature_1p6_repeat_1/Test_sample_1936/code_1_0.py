import itertools

def find_smallest_non_group_n():
    """
    Searches for the smallest integer n for which a non-group structure exists.
    """
    n = 1
    while True:
        elements = range(n)
        found_non_group = False
        
        # A binary operation can be represented as a Cayley table.
        # The number of possible tables is n^(n^2).
        # We can represent the n x n table as a flat tuple of n*n elements.
        num_products = n * n
        
        # Iterate through all possible binary operations for size n.
        # itertools.product generates the Cartesian product, creating all possible tables.
        for flat_table in itertools.product(elements, repeat=num_products):
            # Reconstruct the n x n table from the flat tuple
            table = [list(flat_table[i*n : (i+1)*n]) for i in range(n)]
            
            # Check if this operation forms a group
            is_a_group, reason = check_group(table, n)
            
            if not is_a_group:
                print(f"The smallest number n is {n}.")
                print("An example of a binary operation on a set of {n} elements that is not a group is given by this Cayley table:")
                # Print the Cayley table
                header = "  | " + " ".join(map(str, elements))
                print(header)
                print("---" + "-" * len(header))
                for i, row in enumerate(table):
                    print(f"{i} | " + " ".join(map(str, row)))
                print(f"\nThis structure fails to be a group because: {reason}")
                found_non_group = True
                break
        
        if found_non_group:
            return n # Found the smallest n, so we are done
            
        n += 1

def check_group(table, n):
    """
    Checks if a given Cayley table represents a group.
    
    Returns:
        (bool, str): A tuple containing a boolean indicating if it's a group,
                     and a string with the reason if it's not.
    """
    # Closure is implicitly satisfied by how the table is generated.
    
    # 1. Check for an identity element
    identity_element = None
    for e in range(n):
        is_identity = True
        for i in range(n):
            if table[i][e] != i or table[e][i] != i:
                is_identity = False
                break
        if is_identity:
            identity_element = e
            break
            
    if identity_element is None:
        return False, "No identity element exists."
        
    # 2. Check for an inverse for each element
    for a in range(n):
        has_inverse = False
        for b in range(n):
            # Check if b is the inverse of a
            if table[a][b] == identity_element and table[b][a] == identity_element:
                has_inverse = True
                break
        if not has_inverse:
            return False, f"Element {a} does not have an inverse."
            
    # 3. Check for associativity: (a*b)*c == a*(b*c)
    for a in range(n):
        for b in range(n):
            for c in range(n):
                val1 = table[table[a][b]][c]
                val2 = table[a][table[b][c]]
                if val1 != val2:
                    return False, f"Associativity fails, e.g., for ({a}*{b})*{c} != {a}*({b}*{c})."

    return True, "It is a group."

if __name__ == '__main__':
    find_smallest_non_group_n()
