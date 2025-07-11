import math
from collections import defaultdict

def get_I(k):
    """Calculates the number of indecomposable cool strings for k symbols."""
    if k == 0:
        return 1 # Base case for calculations
    if k == 1:
        return 1
    # For k >= 2, the formula is 2^(k-1) * (k-1)!
    return (2**(k - 1)) * math.factorial(k - 1)

# Memoization cache for partition generation
memo_partitions = {}

def get_partitions(s):
    """Generates all unique set partitions for a given frozenset s."""
    # Using frozenset for memoization keys
    if not s:
        # Base case: the only partition of an empty set is a list containing an empty part
        return [[]] 
    if s in memo_partitions:
        return memo_partitions[s]

    element = next(iter(s))
    rest = s - {element}
    sub_partitions = get_partitions(rest)
    
    result = []
    for part in sub_partitions:
        # Case 1: element gets its own new block (part)
        result.append(part + [{element}])
        # Case 2: element is added to one of the existing blocks
        for i in range(len(part)):
            new_part = part[:]
            # create a new set for the modified block
            new_part[i] = part[i].union({element}) 
            result.append(new_part)
            
    memo_partitions[s] = result
    return result

def solve_cool_strings():
    """
    Calculates the number of cool strings of maximal length 3n for a given n.
    """
    try:
        n_str = input("Enter the number of symbols (n): ")
        n = int(n_str)
        if n < 1:
            print("Error: n must be a positive integer.")
            return
    except ValueError:
        print("Error: Invalid input. Please enter an integer.")
        return

    # Generate all partitions of the set {1, ..., n}
    initial_set = frozenset(range(1, n + 1))
    partitions = get_partitions(initial_set)

    # Group partitions by their shape (e.g., for n=4, shape (2,2) for {{1,2},{3,4}})
    # to structure the output explanation.
    partitions_by_shape = defaultdict(int)
    for p in partitions:
        # The shape is a tuple of the sorted sizes of the blocks.
        shape = tuple(sorted([len(block) for block in p], reverse=True))
        partitions_by_shape[shape] += 1
        
    print(f"\nTo find the number of cool strings for n={n}, we sum the contributions from all partitions of the {n} symbols.")
    print("A partition's contribution depends on its shape (the sizes of its blocks).")
    print("The number of ways to form an indecomposable string on k symbols is I(k).")
    print(f"We use I(1) = {get_I(1)}, and for k >= 2, I(k) = 2^(k-1) * (k-1)!")
    print("\n" + "="*20 + " CALCULATION " + "="*20)

    total_cool_strings = 0
    final_sum_terms = []

    # Iterate through sorted shapes for a consistent output order
    for shape, num_partitions_of_shape in sorted(partitions_by_shape.items()):
        
        num_blocks = len(shape)
        
        # Calculate product of I(k) for each block size k in the shape
        prod_I = 1
        i_parts_str = []
        for k in shape:
            ik_val = get_I(k)
            prod_I *= ik_val
            i_parts_str.append(f"I({k})")

        # Contribution for a single partition of this shape
        term_val = math.factorial(num_blocks) * prod_I
        
        # Total contribution for all partitions of this shape
        total_contribution = num_partitions_of_shape * term_val
        
        final_sum_terms.append(total_contribution)
        total_cool_strings += total_contribution

        # Print the breakdown for this shape
        print(f"\nShape {shape}:")
        print(f"  - Number of partitions with this shape: {num_partitions_of_shape}")
        print(f"  - Contribution per partition = (num_blocks!) * product(I values)")
        print(f"    = {num_blocks}! * ( {' * '.join(i_parts_str)} ) = {term_val}")
        print(f"  - Subtotal for this shape = {num_partitions_of_shape} * {term_val} = {total_contribution}")
    
    print("\n" + "="*20 + " FINAL RESULT " + "="*20)
    equation_str = " + ".join(map(str, final_sum_terms))
    print(f"The total number of cool strings is the sum of these subtotals:")
    if len(final_sum_terms) > 1:
        print(f"Total = {equation_str} = {total_cool_strings}")
    else:
        print(f"Total = {total_cool_strings}")

# Execute the solution
solve_cool_strings()