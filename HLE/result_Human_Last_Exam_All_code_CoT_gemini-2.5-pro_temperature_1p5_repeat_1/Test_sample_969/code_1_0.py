def solve_sequence():
    """
    This function completes the given sequence based on an identified pattern.

    The sequence is grouped into blocks of three:
    B1 = [3, 2, 1]
    B2 = [2, 3, 3]
    B3 = [3, 2, 2]

    The rule to generate the next block B(i+1) from B(i)=[x,y,z] is:
    - The first two elements are swapped: new_x = y, new_y = x.
    - The third element, z, follows a cycling pattern: 1, 3, 2, 1, 3, 2, ...

    Applying the rule:
    - B4 is derived from B3=[3,2,2]:
      - new_x = 2
      - new_y = 3
      - new_z is the next in the cycle, which is 1.
      - So, B4 = [2, 3, 1].

    The next three elements are 2, 3, 1. Option C is [2, 3, 1, 2], which is the closest match.
    We will use the elements from option C to complete the sequence.
    """
    sequence = [3, 2, 1, 2, 3, 3, 3, 2, 2]
    next_four_elements = [2, 3, 1, 2]

    completed_sequence = sequence + next_four_elements

    print("The original sequence is:")
    print(*sequence)
    print("\nThe next 4 elements are:")
    print(*next_four_elements)
    print("\nThe completed sequence is:")
    print(*completed_sequence)

solve_sequence()