import sys

def solve_sequence():
    """
    Finds the next three terms in the sequence based on the rule of
    Pythagorean triples modulo 26.
    
    The rule is: val(L1)^2 + val(L2)^2 is congruent to val(L3)^2 (modulo 26).
    val(A) = 0, val(B) = 1, ..., val(Z) = 25.
    """
    
    # The last letter group in the sequence starts with 'N', so the next is 'O'.
    first_letter_val = ord('O') - ord('A')
    
    # Pre-compute squares modulo 26 for all possible letter values (0-25)
    # to make the main loop more efficient.
    squares_mod_26 = {i: (i * i) % 26 for i in range(26)}
    
    first_letter_sq = squares_mod_26[first_letter_val]
    
    found_terms = []
    
    # Iterate through all possible second letters (from 'A' to 'Z').
    for second_letter_val in range(26):
        second_letter_sq = squares_mod_26[second_letter_val]
        
        # Calculate the left side of the congruence relation.
        left_side = (first_letter_sq + second_letter_sq) % 26
        
        # Iterate through all possible third letters to find a match.
        for third_letter_val in range(26):
            third_letter_sq = squares_mod_26[third_letter_val]
            
            # Check if the rule is satisfied.
            if left_side == third_letter_sq:
                # Construct the three-letter code from the numerical values.
                l1 = chr(ord('A') + first_letter_val)
                l2 = chr(ord('A') + second_letter_val)
                l3 = chr(ord('A') + third_letter_val)
                term = f"{l1}{l2}{l3}"
                found_terms.append(term)

    # Print the first three found terms, as requested.
    for i in range(min(3, len(found_terms))):
        term = found_terms[i]
        n1 = ord(term[0]) - ord('A')
        n2 = ord(term[1]) - ord('A')
        n3 = ord(term[2]) - ord('A')
        
        # Output the term and the equation with its numbers, as per the instruction.
        print(f"Result: {term}, based on the equation: {n1}^2 + {n2}^2 = {n3}^2 (mod 26)")

solve_sequence()
