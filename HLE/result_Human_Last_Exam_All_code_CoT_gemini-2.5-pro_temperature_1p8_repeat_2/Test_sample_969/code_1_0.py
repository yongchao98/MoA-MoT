def generate_sequence():
    """
    This function generates the next 4 elements of the sequence based on a discovered pattern.
    
    The pattern is identified as a series of groups:
    - ud(N): A 'down-then-up' group, counting from N down to 1, then from 2 up to 3.
    - rep(N): A 'repeat' group, where the number N is repeated N times.
    
    The given sequence `3 2 1 2 3 3 3 2 2` corresponds to:
    - ud(3): [3, 2, 1, 2, 3]
    - rep(3): [3, 3, 3]
    - rep(2): [2, 2]
    
    The pattern implies that after the groups for N=3 are generated, a new cycle begins for N-1=2.
    The next group to be generated is therefore ud(2).
    """
    
    # Generate the parts for the ud(2) group
    # Part 1: Count down from 2 to 1
    down_part = [2, 1]
    
    # Part 2: Count up from 2 to 3
    up_part = [2, 3]
    
    # The next four elements are the concatenation of these two parts
    next_elements = down_part + up_part
    
    print("The given sequence is: 3 2 1 2 3 3 3 2 2")
    print("The pattern continues by generating a 'down-then-up' group starting from 2.")
    print(f"Down from 2 to 1 is: {down_part[0]} {down_part[1]}")
    print(f"Up from 2 to 3 is: {up_part[0]} {up_part[1]}")
    print("Combining them gives the next 4 elements:")
    
    # Print each number of the result, as requested
    print(f"{next_elements[0]} {next_elements[1]} {next_elements[2]} {next_elements[3]}")

generate_sequence()