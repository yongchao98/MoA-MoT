def find_next_codes():
    """
    Calculates the next three codes in the sequence based on the rule
    that the letter positions form an arithmetic progression.
    """
    # The last letter in the sequence is 'N', so the next is 'O'.
    x_char = 'O'
    x = ord(x_char) - ord('A')
    
    # We will find the codes for the first three possible second letters: 'A', 'B', 'C'.
    y_chars = ['A', 'B', 'C']
    
    next_codes = []
    
    print("The rule is: Z = (2*Y - X) mod 26, where X, Y, Z are 0-indexed letter positions.")
    print("-" * 20)
    
    for y_char in y_chars:
        y = ord(y_char) - ord('A')
        
        # Apply the arithmetic progression formula
        z = (2 * y - x) % 26
        
        z_char = chr(ord('A') + z)
        code = f"{x_char}{y_char}{z_char}"
        next_codes.append(code)
        
        print(f"For X='{x_char}' (x={x}) and Y='{y_char}' (y={y}):")
        print(f"  z = (2 * {y} - {x}) mod 26 = {z}")
        print(f"  The resulting code is: {code}")
        print("-" * 20)

    print(f"\nThe next three capital letters in the sequence are: {', '.join(next_codes)}")

find_next_codes()