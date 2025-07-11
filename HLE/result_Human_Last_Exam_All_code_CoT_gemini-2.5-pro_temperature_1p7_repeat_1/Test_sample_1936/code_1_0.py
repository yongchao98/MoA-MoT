import itertools

def is_group(n, table):
    """
    Checks if a given Cayley table of size n represents a group.
    The set is assumed to be {0, 1, ..., n-1}.
    """
    # 1. Check Associativity
    for i in range(n):
        for j in range(n):
            for k in range(n):
                if table[table[i][j]][k] != table[i][table[j][k]]:
                    return False  # Associativity failed

    # 2. Find Identity Element
    identity_element = None
    for e in range(n):
        is_identity = True
        for a in range(n):
            if table[e][a] != a or table[a][e] != a:
                is_identity = False
                break
        if is_identity:
            identity_element = e
            break
    
    if identity_element is None:
        return False # No identity element found

    # 3. Check for Inverse Element for each element
    for a in range(n):
        has_inverse = False
        for b in range(n):
            if table[a][b] == identity_element and table[b][a] == identity_element:
                has_inverse = True
                break
        if not has_inverse:
            return False # Element 'a' has no inverse

    return True

def find_smallest_n_for_non_group():
    """
    Finds the smallest positive integer n for which a set of size n
    can have a binary operation that does not form a group.
    """
    n = 1
    while True:
        elements = range(n)
        found_non_group = False

        # Generate all possible binary operations (Cayley tables)
        # A table is an n x n matrix where each cell can be any of n elements.
        # Total number of tables is n^(n*n).
        
        # Generator for all possible rows in the Cayley table
        possible_rows = itertools.product(elements, repeat=n)
        # Generator for all possible tables (combinations of rows)
        all_tables = itertools.product(possible_rows, repeat=n)

        for table_tuple in all_tables:
            table = [list(row) for row in table_tuple]
            if not is_group(n, table):
                print(f"The smallest positive integer n for which a non-group exists is:")
                # The question asks to "output each number in the final equation!".
                # We will interpret this as printing the final answer clearly.
                print(n)
                
                print("\nAn example of a binary operation on a set of this size that is not a group is given by the following Cayley table:")
                header = "  . | " + " ".join(map(str, elements))
                print(header)
                print("---" + "+---" * n)
                for i, row in enumerate(table):
                    row_str = f"  {i} | " + " ".join(map(str, row))
                    print(row_str)

                found_non_group = True
                break
        
        if found_non_group:
            break
        
        n += 1

if __name__ == '__main__':
    find_smallest_n_for_non_group()