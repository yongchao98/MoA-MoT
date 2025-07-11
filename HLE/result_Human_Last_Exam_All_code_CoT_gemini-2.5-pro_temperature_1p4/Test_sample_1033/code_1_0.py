def solve_sequence():
    """
    Calculates the next three terms in the sequence based on a plausible pattern.
    """
    # The first letter of the next three terms.
    l1_char = 'O'
    
    # Assume the second letter starts at 'A' and increments.
    next_l2_chars = ['A', 'B', 'C']
    
    print("Based on the likely intended pattern (L1 + L2 = L3), the next three terms are:")
    
    for l2_char in next_l2_chars:
        # Convert letters to numbers 1-26
        p1 = ord(l1_char) - ord('A') + 1
        p2 = ord(l2_char) - ord('A') + 1
        
        # Alphabetic addition: (p1 + p2 - 1) mod 26 + 1
        # The -1 and +1 handle the wrap-around for a 1-26 range.
        p3 = (p1 + p2 - 1) % 26 + 1
        
        # Convert result back to a character
        l3_char = chr(p3 + ord('A') - 1)
        
        # Output the equation and the resulting triplet
        print(f"{l1_char} ({p1}) + {l2_char} ({p2}) = {l3_char} ({p3})  --> Triplet: {l1_char}{l2_char}{l3_char}")

solve_sequence()