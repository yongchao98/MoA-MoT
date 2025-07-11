def solve_sequence():
    """
    This function calculates the next three terms in the given cryptic sequence.
    """
    # The last group in the sequence is 'N', so the next group is 'O'.
    l1 = 'O'
    v1 = ord(l1) - ord('A')

    # The second letters for the 'O' group are predicted to be 'C', 'E', 'H'
    next_l2s = ['C', 'E', 'H']
    
    # The formula deduced is: v3 = (v1 * 5 + v2 * 17 + i * 17 + 24) % 26
    # where i is the index in the group.
    
    print("Calculating the next three terms for the sequence:")
    
    for i, l2 in enumerate(next_l2s):
        v2 = ord(l2) - ord('A')
        
        # Apply the formula
        a = v1 * 5
        b = v2 * 17
        c = i * 17
        d = 24
        v3 = (a + b + c + d) % 26
        l3 = chr(v3 + ord('A'))
        
        # Print the equation and the result
        print(f"Term {i+1}:")
        print(f"L1='{l1}'(v1={v1}), L2='{l2}'(v2={v2}), i={i}")
        print(f"v3 = ({v1} * 5 + {v2} * 17 + {i} * 17 + 24) % 26")
        print(f"v3 = ({a} + {b} + {c} + {d}) % 26")
        v3_calc = (a+b+c+d)
        print(f"v3 = {v3_calc} % 26 = {v3}")
        print(f"L3='{l3}'")
        print(f"Result: {l1}{l2}{l3}\n")

solve_sequence()
