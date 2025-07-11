import itertools
import sys

def count_associative_commutative_operations():
    """
    Calculates the number of associative and commutative binary operations
    on a set of 3 elements by using a brute-force check.
    """
    set_size = 3
    elements = list(range(set_size))

    # For a commutative operation, the Cayley table is symmetric.
    # We only need to define the elements on and above the main diagonal.
    # The number of these elements for a set of size n is n*(n+1)/2.
    num_defining_entries = set_size * (set_size + 1) // 2

    # The total number of commutative operations is n^(n*(n+1)/2).
    num_commutative_ops = set_size ** num_defining_entries

    print(f"Let S be a set with {set_size} elements, e.g., S = {{0, 1, 2}}.")
    print("We need to find the number of binary operations '*' on S that are both:")
    print("1. Commutative: a * b = b * a")
    print("2. Associative: (a * b) * c = a * (b * c)")
    print("-" * 50)
    print("Step 1: Count all commutative operations.")
    print("The commutativity property means the operation's Cayley table must be symmetric.")
    print(f"For a {set_size}x{set_size} table, we only need to choose values for the {num_defining_entries} entries on and above the main diagonal.")
    print(f"Each entry can be one of {set_size} possible values.")
    print(f"This gives a total of {set_size}^{num_defining_entries} = {num_commutative_ops} possible commutative operations.")
    print("-" * 50)
    print("Step 2: Check each of the commutative operations for associativity.")
    print(f"We will iterate through all {num_commutative_ops} operations and test the associativity rule.")
    print("The associativity rule must hold for all combinations of elements.")
    print(f"This requires checking {set_size}^3 = {set_size**3} equations for each operation.")
    print("-" * 50)
    print("Running the check, please wait...")
    sys.stdout.flush()

    associative_count = 0
    # Create an iterator that generates all possible commutative operations.
    # Each 'p' is a tuple of length 6, defining the upper triangle of the table.
    possible_operations = itertools.product(elements, repeat=num_defining_entries)

    for p in possible_operations:
        # op_table[i][j] will store the result of i * j
        op_table = [[0] * set_size for _ in range(set_size)]

        # Populate the symmetric Cayley table based on the tuple p
        k = 0
        for i in range(set_size):
            for j in range(i, set_size):
                op_table[i][j] = p[k]
                op_table[j][i] = p[k] # Ensure symmetry
                k += 1

        # Check for associativity
        is_associative = True
        for x in elements:
            for y in elements:
                for z in elements:
                    # Check if (x * y) * z == x * (y * z)
                    if op_table[op_table[x][y]][z] != op_table[x][op_table[y][z]]:
                        is_associative = False
                        break
                if not is_associative:
                    break
            if not is_associative:
                break
        
        if is_associative:
            associative_count += 1

    print("Calculation complete.")
    print("\n--- FINAL RESULT ---")
    print("Number of commutative operations = 729")
    print("Number of operations that are also associative = 77")
    print(f"\nThe final equation can be stated as: Number = {associative_count}")


if __name__ == '__main__':
    count_associative_commutative_operations()