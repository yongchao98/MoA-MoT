import itertools

def is_group(n, table):
    """
    Checks if a given Cayley table of size n represents a group.
    Returns a tuple (is_group, reason_for_failure).
    """
    # Closure is guaranteed by how the table is generated.

    # 1. Check Associativity: (a.b).c == a.(b.c)
    for a in range(n):
        for b in range(n):
            for c in range(n):
                if table[table[a][b]][c] != table[a][table[b][c]]:
                    return False, "Failed Associativity"

    # 2. Find Identity Element e: a.e == e.a == a
    identity = None
    for e in range(n):
        is_identity = True
        for a in range(n):
            if table[a][e] != a or table[e][a] != a:
                is_identity = False
                break
        if is_identity:
            identity = e
            break
    
    if identity is None:
        return False, "Failed to find an Identity Element"

    # 3. Check for Inverse Element b for each a: a.b == b.a == e
    for a in range(n):
        has_inverse = False
        for b in range(n):
            if table[a][b] == identity and table[b][a] == identity:
                has_inverse = True
                break
        if not has_inverse:
            return False, f"Failed to find an Inverse for element {a}"
            
    return True, "Forms a group"

def find_smallest_n_for_non_group():
    """
    Finds the smallest n >= 1 for which a non-group structure exists.
    """
    n = 1
    while True:
        is_always_group = True
        
        num_entries = n * n
        elements = range(n)
        
        # Iterate through all possible binary operations (Cayley tables) for size n
        for flat_table in itertools.product(elements, repeat=num_entries):
            table = [list(flat_table[i*n : (i+1)*n]) for i in range(n)]
            
            is_g, reason = is_group(n, table)
            
            if not is_g:
                # Found a non-group. This n is the answer.
                print(f"The smallest number n for which a non-group structure exists is: {n}")
                print("\nAn example of such a structure is on the set G = {" + ", ".join(map(str, range(n))) + "}")
                print("with the binary operation defined by the following Cayley table:")
                
                # Print table header
                header = "  . | " + " | ".join(map(str, range(n)))
                print(header)
                print("----+" + "---" * n)
                
                # Print table rows
                for i in range(n):
                    row_str = f"  {i} | " + " | ".join(map(str, table[i]))
                    print(row_str)
                
                print(f"\nThis structure is not a group. Reason: {reason}.")
                
                print("\nThe defining equations for this operation are:")
                for i in range(n):
                    for j in range(n):
                        # This satisfies "output each number in the final equation!"
                        print(f"{i} . {j} = {table[i][j]}")
                
                is_always_group = False
                return # Exit after finding the first n

        # If we reach here, it means all operations for size n formed a group.
        # This will happen for n=1.
        n += 1

if __name__ == '__main__':
    find_smallest_n_for_non_group()