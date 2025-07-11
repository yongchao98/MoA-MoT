def solve_sequence():
    """
    This function calculates the next three 3-letter codes in the sequence
    based on a set of polynomial formulas.
    """
    
    # The first letter of the next group is 'O'.
    v1 = 14  # O is the 15th letter, so its 0-index is 14.
    
    # Formulas found to generate the sequence
    # v2 = (14*v1^2 + 5*k^2 + 7*v1*k + 22*v1 + 11*k + 2) mod 26
    # v3 = (11*v1 + 17*v2 + 13*k + 19) mod 26
    # where v1 is the value of the first letter (A=0), and
    # k is the index within the group.

    results = []
    print("The next three terms are derived from the following equations:")
    print("v2 = (14*v1^2 + 5*k^2 + 7*v1*k + 22*v1 + 11*k + 2) % 26")
    print("v3 = (11*v1 + 17*v2 + 13*k + 19) % 26\n")

    # We need to find the next 3 terms, so k will be 0, 1, 2
    for k in range(3):
        # Coefficients from the formula for v2
        c1_v2, c2_v2, c3_v2, c4_v2, c5_v2, c6_v2 = 14, 5, 7, 22, 11, 2

        # Calculate v2
        v2 = (c1_v2*v1**2 + c2_v2*k**2 + c3_v2*v1*k + c4_v2*v1 + c5_v2*k + c6_v2) % 26
        
        # Coefficients from the formula for v3
        c1_v3, c2_v3, c3_v3, c4_v3 = 11, 17, 13, 19
        
        # Calculate v3
        v3 = (c1_v3*v1 + c2_v3*v2 + c3_v3*k + c4_v3) % 26
        
        # Convert v1, v2, v3 back to letters
        l1 = chr(ord('A') + v1)
        l2 = chr(ord('A') + v2)
        l3 = chr(ord('A') + v3)
        
        # Print the derivation for each term
        print(f"For the term with k={k}:")
        print(f"v1 = {v1}")
        print(f"v2 = ({c1_v2}*{v1}^2 + {c2_v2}*{k}^2 + {c3_v2}*{v1}*{k} + {c4_v2}*{v1} + {c5_v2}*{k} + {c6_v2}) % 26 = {v2}")
        print(f"v3 = ({c1_v3}*{v1} + {c2_v3}*{v2} + {c3_v3}*{k} + {c4_v3}) % 26 = {v3}")
        print(f"Resulting term: {l1}{l2}{l3}\n")
        
        results.append(f"{l1}{l2}{l3}")
        
    final_answer = " ".join(results)
    # The final answer format is specific
    print(f"<<<{final_answer}>>>")

solve_sequence()