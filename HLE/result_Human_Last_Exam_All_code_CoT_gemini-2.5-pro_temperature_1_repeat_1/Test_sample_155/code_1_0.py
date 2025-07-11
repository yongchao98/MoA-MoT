import collections

def solve_braille_puzzle():
    """
    Decodes the word from the puzzle based on the Braille encoding scheme.
    """
    # Step 1: Parse the input strings from the puzzle
    top_str = "t b b t t t t b b b"
    mid_str = "t b t t b b b b b b"
    bot_str = "b b t t t t b t b b"

    # Convert strings to lists of 't' and 'b'
    top_h = top_str.split()
    mid_h = mid_str.split()
    bot_h = bot_str.split()

    # The length of the horizontal strings determines the number of characters to decode.
    # In this case, 10 characters, but the actual word might be shorter.
    # The symmetrical decoding suggests a word of length 5 (10/2). Let's test this.
    word_len = len(top_h) // 2

    # Braille alphabet mapping (6-bit binary to character)
    # Bits are ordered: 1 2 3 4 5 6
    braille_map = {
        "100000": "a", "110000": "b", "100100": "c", "100110": "d",
        "100010": "e", "110100": "f", "110110": "g", "110010": "h",
        "010100": "i", "010110": "j", "101000": "k", "111000": "l",
        "101100": "m", "101110": "n", "101010": "o", "111100": "p",
        "111110": "q", "111010": "r", "011100": "s", "011110": "t",
        "101001": "u", "111001": "v", "010111": "w", "101101": "x",
        "101111": "y", "101011": "z", "000000": " "
    }

    decoded_word = []
    
    # Final Equation will show the decoding for the first letter
    final_equation = "First letter decoding:\n"
    
    # Step 2 & 3: Decode each letter
    for i in range(word_len):
        # Symmetrical mapping: char i is paired with char (length - 1 - i)
        j = len(top_h) - 1 - i

        # Map 't' to '1' (dot) and 'b' to '0' (no dot)
        # Left column of dots (1, 2, 3) from left-to-right
        dot1 = '1' if top_h[i] == 't' else '0'
        dot2 = '1' if mid_h[i] == 't' else '0'
        dot3 = '1' if bot_h[i] == 't' else '0'

        # Right column of dots (4, 5, 6) from right-to-left
        dot4 = '1' if top_h[j] == 't' else '0'
        dot5 = '1' if mid_h[j] == 't' else '0'
        dot6 = '1' if bot_h[j] == 't' else '0'

        braille_code = dot1 + dot2 + dot3 + dot4 + dot5 + dot6
        
        # Store the equation for the first letter
        if i == 0:
            final_equation += f"Dot 1 from top[{i}]: {dot1}\n"
            final_equation += f"Dot 2 from mid[{i}]: {dot2}\n"
            final_equation += f"Dot 3 from bot[{i}]: {dot3}\n"
            final_equation += f"Dot 4 from top[{j}]: {dot4}\n"
            final_equation += f"Dot 5 from mid[{j}]: {dot5}\n"
            final_equation += f"Dot 6 from bot[{j}]: {dot6}\n"
            final_equation += f"Braille code: {dot1}{dot2}{dot3}{dot4}{dot5}{dot6}\n"
            char = braille_map.get(braille_code, '?').upper()
            final_equation += f"Character: '{char}'\n"

        # Look up character in map, convert to upper case for the final word
        decoded_char = braille_map.get(braille_code, '?')
        decoded_word.append(decoded_char.upper())

    print(final_equation)
    print("The decoded word is:")
    print("".join(decoded_word))

solve_braille_puzzle()