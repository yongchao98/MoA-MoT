def solve_sequence():
    """
    Solves the sequence puzzle by identifying and extending the pattern.
    """
    initial_sequence = [3, 2, 1, 2, 3, 3, 3, 2, 2]
    
    # The pattern is determined to be:
    # 1. A sweep from 3 down to 1 and back to 3: 3, 2, 1, 2, 3
    # 2. A block of three 3s: 3, 3, 3
    # 3. A block of two 2s: 2, 2
    # This concludes "Stage 3". The next stage, "Stage 2", begins.
    # 4. The next part is a sweep starting from 2, going down to 1, and then up to the previous stage number (3): 2, 1, 2, 3
    
    next_elements = [2, 1, 2, 3]
    
    full_sequence = initial_sequence + next_elements
    
    print("The initial sequence is: ", end="")
    print(*initial_sequence, sep=' ')
    
    print("The pattern continues with a new sweep (2 -> 1 -> 3).")
    
    print("The next 4 elements are: ", end="")
    print(*next_elements, sep=' ')
    
    print("The completed sequence is: ", end="")
    print(*full_sequence, sep=' ')

solve_sequence()