def solve_sequence():
    """
    This function generates the sequence based on the derived block pattern
    and prints the next 4 elements.
    """
    # The given prefix of the sequence
    given_sequence = [3, 2, 1, 2, 3, 3, 3, 2, 2]

    # Pattern for the 'values' sequence (v_k)
    # It bounces between 3 and 1.
    v = []
    current_v = 3
    direction = -1
    for _ in range(10): # Generate enough terms
        v.append(current_v)
        if current_v == 1:
            direction = 1
        elif current_v == 3:
            direction = -1
        current_v += direction
        
    # Pattern for the 'counts' sequence (n_k)
    # It repeats the pattern (1, 1, 1, 1, 3, 3)
    n_pattern = [1, 1, 1, 1, 3, 3]
    n = []
    for i in range(10): # Generate enough terms
        n.append(n_pattern[i % len(n_pattern)])

    # Generate the full sequence from the block patterns
    full_sequence = []
    blocks = []
    for i in range(len(v)):
        block = [str(v[i])] * n[i]
        blocks.append(" ".join(block))
        for j in range(n[i]):
            full_sequence.append(v[i])
            
    # The next 4 elements are the ones after the given 9-element sequence
    next_elements = full_sequence[len(given_sequence):len(given_sequence) + 4]

    print("The given sequence is a prefix of a sequence constructed from blocks.")
    print("Underlying value sequence (v): {}".format(v[:9]))
    print("Underlying count sequence (n): {}".format(n[:9]))
    print("\nThis generates the following blocks:")
    print(" -> ".join([f"({s})" for s in blocks[:9]]))

    print("\nOriginal sequence: {}".format(" ".join(map(str, given_sequence))))
    print("The next 4 elements are: {}".format(" ".join(map(str, next_elements))))
    
    print("\nThe full sequence starts: {} ...".format(" ".join(map(str, full_sequence[:len(given_sequence)+4]))))
    
    # Final equation format as requested
    print("\nFinal Equation:")
    print("{} + {} = {}".format(
        " ".join(map(str, given_sequence)),
        " ".join(map(str, next_elements)),
        " ".join(map(str, full_sequence[:len(given_sequence)+4]))
    ))

solve_sequence()