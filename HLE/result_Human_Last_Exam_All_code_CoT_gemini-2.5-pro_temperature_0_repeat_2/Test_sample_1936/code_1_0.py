import itertools

def find_smallest_non_group():
    """
    Finds the smallest integer n >= 1 for which there exists a binary operation
    on a set of n elements that does not form a group.
    """

    # Start checking from n=1, as groups are conventionally non-empty.
    n = 1
    while True:
        print(f"--- Checking for n = {n} ---")
        elements = list(range(n))
        num_operations = n**(n*n)
        
        is_always_a_group = True
        
        # Generate all possible binary operations (Cayley tables) for size n.
        # An operation is a function GxG -> G, which can be represented by an n x n table.
        # There are n choices for each of the n*n entries in the table.
        # We use itertools.product to generate all possible tables.
        for table_flat in itertools.product(elements, repeat=n*n):
            table = [list(table_flat[i*n : (i+1)*n]) for i in range(n)]
            
            # Check if this operation table forms a group.
            is_g, reason, equation = check_group_properties(table, n)
            
            if not is_g:
                # We found an operation for this n that is NOT a group.
                # This means n is our answer.
                is_always_a_group = False
                print(f"\nSuccess! Found a non-group structure for n = {n}.")
                print("This is the smallest size for a non-empty set with a binary operation that is not a group.")
                print("\nCounterexample Cayley Table (elements are {0, 1, ...}):")
                header = "  . | " + " ".join(map(str, elements))
                print(header)
                print("----+" + "-" * len(header))
                for i, row in enumerate(table):
                    print(f"  {i} | " + " ".join(map(str, row)))
                
                print(f"\nReason for failure: {reason}")
                print(f"Verification: {equation}")
                
                print(f"\nTherefore, the smallest number n is {n}.")
                return n

        if is_always_a_group:
            print(f"All {num_operations} binary operations for n = {n} form a group.")
        
        # If all operations for the current n form a group, try the next n.
        n += 1

def check_group_properties(table, n):
    """
    Checks if a given Cayley table for a set of n elements forms a group.
    Returns (is_group, reason, equation).
    """
    elements = list(range(n))

    # 1. Check for an identity element
    identity = None
    for e_cand in elements:
        is_identity = True
        for x in elements:
            if table[e_cand][x] != x or table[x][e_cand] != x:
                is_identity = False
                break
        if is_identity:
            identity = e_cand
            break
    
    if identity is None:
        # Provide a specific counterexample for why no element is the identity
        e_cand = elements[0] # Check the first candidate
        for x in elements:
            if table[e_cand][x] != x:
                return (False, "No identity element", f"For candidate e={e_cand}, the equation e . x = x fails for x={x}. We have {e_cand} . {x} = {table[e_cand][x]}, which is not {x}.")
            if table[x][e_cand] != x:
                 return (False, "No identity element", f"For candidate e={e_cand}, the equation x . e = x fails for x={x}. We have {x} . {e_cand} = {table[x][e_cand]}, which is not {x}.")
        return (False, "No identity element", "No element satisfies the identity property.")

    # 2. Check for associativity: (a.b).c == a.(b.c)
    for a in elements:
        for b in elements:
            for c in elements:
                ab = table[a][b]
                bc = table[b][c]
                lhs = table[ab][c]
                rhs = table[a][bc]
                if lhs != rhs:
                    equation = f"({a} . {b}) . {c} = {ab} . {c} = {lhs}, but {a} . ({b} . {c}) = {a} . {bc} = {rhs}"
                    return (False, "Fails associativity", equation)

    # 3. Check for an inverse element for all elements
    for a in elements:
        has_inverse = False
        for b in elements:
            if table[a][b] == identity and table[b][a] == identity:
                has_inverse = True
                break
        if not has_inverse:
            equation = f"For element a={a}, there is no b such that a . b = b . a = {identity} (the identity)."
            return (False, f"No inverse for element {a}", equation)
            
    return (True, "Is a group", "")

if __name__ == '__main__':
    find_smallest_non_group()