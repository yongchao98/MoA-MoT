def solve_sequence():
    """
    Calculates the next three terms in the sequence based on a derived formula.
    """
    # Helper to convert a character to its 0-25 value
    def c_to_i(c):
        return ord(c) - ord('A')

    # Helper to convert a 0-25 value back to a character
    def i_to_c(i):
        return chr(i + ord('A'))

    # The derived formula for the third letter's value (v3)
    # v3 = (10*v1 + 2*v2 + 18*j + 7) % 26
    # v1: value of the first letter (L1)
    # v2: value of the second letter (L2)
    # j: 0-based index of the term within its L1-block
    def calculate_v3(v1, v2, j):
        return (10 * v1 + 2 * v2 + 18 * j + 7) % 26

    # The next block is 'O'
    v1 = c_to_i('O')
    
    # We assume the sequence of second letters (L2) for 'O' repeats the 'A' block's sequence
    l2_sequence_for_a = ['C', 'C', 'D']
    
    next_terms = []
    
    # Calculate the first term (j=0)
    j = 0
    l1 = i_to_c(v1)
    l2 = l2_sequence_for_a[j]
    v2 = c_to_i(l2)
    v3 = calculate_v3(v1, v2, j)
    l3 = i_to_c(v3)
    next_terms.append(f"{l1}{l2}{l3}")
    
    # Calculate the second term (j=1)
    j = 1
    l1 = i_to_c(v1)
    l2 = l2_sequence_for_a[j]
    v2 = c_to_i(l2)
    v3 = calculate_v3(v1, v2, j)
    l3 = i_to_c(v3)
    next_terms.append(f"{l1}{l2}{l3}")

    # Calculate the third term (j=2)
    j = 2
    l1 = i_to_c(v1)
    l2 = l2_sequence_for_a[j]
    v2 = c_to_i(l2)
    v3 = calculate_v3(v1, v2, j)
    l3 = i_to_c(v3)
    next_terms.append(f"{l1}{l2}{l3}")
    
    print(f"The formula is: L3 = (10 * L1 + 2 * L2 + 18 * j + 7) % 26")
    print(f"Assuming the L2 sequence for 'O' starts like 'A' (C, C, D...).")
    print("-" * 20)
    print("Calculation for the first term:")
    print(f"L1 = 'O'({c_to_i('O')}), L2 = 'C'({c_to_i('C')}), j = 0")
    print(f"L3 = (10 * {c_to_i('O')} + 2 * {c_to_i('C')} + 18 * 0 + 7) % 26 = {i_to_c(calculate_v3(c_to_i('O'), c_to_i('C'), 0))}")
    print("-" * 20)
    print("Calculation for the second term:")
    print(f"L1 = 'O'({c_to_i('O')}), L2 = 'C'({c_to_i('C')}), j = 1")
    print(f"L3 = (10 * {c_to_i('O')} + 2 * {c_to_i('C')} + 18 * 1 + 7) % 26 = {i_to_c(calculate_v3(c_to_i('O'), c_to_i('C'), 1))}")
    print("-" * 20)
    print("Calculation for the third term:")
    print(f"L1 = 'O'({c_to_i('O')}), L2 = 'D'({c_to_i('D')}), j = 2")
    print(f"L3 = (10 * {c_to_i('O')} + 2 * {c_to_i('D')} + 18 * 2 + 7) % 26 = {i_to_c(calculate_v3(c_to_i('O'), c_to_i('D'), 2))}")
    print("-" * 20)
    print(f"The next three capital letters are: {', '.join(next_terms)}")

solve_sequence()