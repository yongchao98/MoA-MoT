import string

def solve_sequence():
    """
    Solves the letter sequence puzzle by finding the underlying mathematical rule
    and generating the next three terms.
    """
    key_map = {
        'A': 2, 'B': 2, 'C': 2, 'D': 3, 'E': 3, 'F': 3, 'G': 4, 'H': 4, 'I': 4,
        'J': 5, 'K': 5, 'L': 5, 'M': 6, 'N': 6, 'O': 6, 'P': 7, 'Q': 7, 'R': 7,
        'S': 7, 'T': 8, 'U': 8, 'V': 8, 'W': 9, 'X': 9, 'Y': 9, 'Z': 9
    }
    
    # Using 0-indexed positions
    letters_on_key = {
        2: ['A', 'B', 'C'], 3: ['D', 'E', 'F'], 4: ['G', 'H', 'I'],
        5: ['J', 'K', 'L'], 6: ['M', 'N', 'O'], 7: ['P', 'Q', 'R', 'S'],
        8: ['T', 'U', 'V'], 9: ['W', 'X', 'Y', 'Z']
    }
    
    pos_map = {
        letter: pos for key, letters in letters_on_key.items() 
        for pos, letter in enumerate(letters)
    }
    
    key_size_map = {key: len(letters) for key, letters in letters_on_key.items()}

    def is_valid(l1, l2, l3):
        """Checks if a triplet (l1, l2, l3) is valid based on the phone keypad rule."""
        k1, k2, k3 = key_map[l1], key_map[l2], key_map[l3]
        p1, p2, p3 = pos_map[l1], pos_map[l2], pos_map[l3]
        
        # Rule for the key
        k3_pred = ((k1 - 2) * (k2 - 2)) % 8 + 2
        
        if k3 != k3_pred:
            return False, None
            
        # Rule for the position
        size_k3 = key_size_map[k3]
        p3_pred = (p1 + p2) % size_k3
        
        if p3 == p3_pred:
            # Prepare the equation string for printing
            eq = (
                f"For {l1}{l2}{l3}:\n"
                f"  Key rule: k({l3}) = ((k({l1})-2) * (k({l2})-2)) % 8 + 2  "
                f"->  {k3} = (({k1}-2) * ({k2}-2)) % 8 + 2 = {k3_pred}\n"
                f"  Position rule: p({l3}) = (p({l1}) + p({l2})) % size(k({l3}))  "
                f"->  {p3} = ({p1} + {p2}) % {size_k3} = {p3_pred}"
            )
            return True, eq
        return False, None

    found_triplets = []
    
    # Start searching alphabetically after 'NZX'
    start_char_code = ord('N')
    
    for i in range(start_char_code, ord('Z') + 1):
        l1 = chr(i)
        for l2 in string.ascii_uppercase:
            # Optimization: determine the required key for l3
            k1, k2 = key_map[l1], key_map[l2]
            k3_target = ((k1 - 2) * (k2 - 2)) % 8 + 2
            
            # Optimization: determine the required position for l3
            p1, p2 = pos_map[l1], pos_map[l2]
            p3_target = (p1 + p2) % key_size_map[k3_target]
            
            l3 = letters_on_key[k3_target][p3_target]
            
            # Ensure we are only looking for triplets after NZX
            triplet_str = f"{l1}{l2}{l3}"
            if triplet_str <= "NZX":
                continue

            is_valid_res, eq_str = is_valid(l1, l2, l3)
            if is_valid_res:
                found_triplets.append((triplet_str, eq_str))
                if len(found_triplets) == 3:
                    break
        if len(found_triplets) == 3:
            break
            
    print("Found the next three triplets in the sequence:\n")
    for triplet, equation in found_triplets:
        print(f"Triplet: {triplet}")
        print(equation)
        print("-" * 20)
        
solve_sequence()