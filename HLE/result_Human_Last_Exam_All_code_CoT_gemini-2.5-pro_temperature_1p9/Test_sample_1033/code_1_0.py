import string

def solve_sequence():
    """
    Solves for the next three terms in the given letter sequence.
    """
    
    # The rule is derived from a linear congruence analysis of the sequence terms.
    # Let v1, v2, v3 be the 0-indexed values of the letters (A=0, ..., Z=25).
    # The equation governing the sequence is:
    # 0*v1 + 0*v2 + 13*v3 = 13 (mod 26)
    # This simplifies to 13*v3 = 13 (mod 26), which holds true if and only if
    # v3 is an odd number.
    
    equation_v1_coeff = 0
    equation_v2_coeff = 0
    equation_v3_coeff = 13
    equation_const = 13
    
    print(f"The rule discovered is a mathematical equation relating the 0-indexed values of the letters (v1, v2, v3):")
    print(f"{equation_v1_coeff}*v1 + {equation_v2_coeff}*v2 + {equation_v3_coeff}*v3 = {equation_const} (mod 26)\n")

    def check_rule(triplet_str):
        v3 = ord(triplet_str[2]) - ord('A')
        return v3 % 2 == 1

    last_term = "NZX"
    
    # The next term must start with the letter after 'N', which is 'O'.
    next_l1_char = chr(ord(last_term[0]) + 1)
    
    found_terms = []
    
    # The sequence is sorted lexicographically. We find the first 3 valid terms
    # that come after the given sequence.
    
    start_l1_val = ord(next_l1_char) - ord('A')
    
    # Iterate from the next main letter (O) onwards
    for v1 in range(start_l1_val, 26):
        if len(found_terms) >= 3:
            break
        l1 = chr(v1 + ord('A'))
        
        # Iterate through the second letter from A to Z
        for l2 in string.ascii_uppercase:
            if len(found_terms) >= 3:
                break
            
            # Iterate through the third letter from A to Z
            for l3 in string.ascii_uppercase:
                triplet = f"{l1}{l2}{l3}"
                v3 = ord(l3) - ord('A')
                # Check if the triplet satisfies the rule
                if (equation_v3_coeff * v3) % 26 == equation_const % 26:
                    found_terms.append(triplet)
                    if len(found_terms) >= 3:
                        break

    print("The next three capital letters in the sequence are:")
    for term in found_terms:
        print(term)

solve_sequence()