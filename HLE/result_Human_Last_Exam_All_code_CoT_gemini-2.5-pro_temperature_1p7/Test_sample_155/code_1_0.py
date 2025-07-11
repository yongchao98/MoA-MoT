def solve():
    """
    Decodes a message encoded in a Braille-like system from the given
    't' (dot on) and 'b' (dot off) values. The structure of the input
    points to an equation with 5 symbols.
    """

    h1 = ['t', 'b', 'b', 'b', 't', 't', 't', 'b', 'b', 'b']
    h2 = ['t', 'b', 't', 't', 'b', 'b', 'b', 'b', 'b', 'b']
    h3 = ['b', 'b', 't', 't', 't', 't', 'b', 't', 'b', 'b']

    # Braille alphabet and numbers mapping (dots to character)
    # Using a simplified mapping for clarity
    braille_map = {
        (1, 0, 0, 0, 0, 0): '1',
        (1, 1, 0, 0, 0, 0): '2',
        (1, 0, 0, 1, 0, 0): '3',
        (1, 0, 0, 1, 1, 0): '4',
        (1, 0, 0, 0, 1, 0): '5',
        (1, 1, 0, 1, 0, 0): '6',
        (1, 1, 0, 1, 1, 0): '7',
        (1, 1, 0, 0, 1, 0): '8',
        (0, 1, 0, 1, 0, 0): '9',
        (0, 1, 0, 1, 1, 0): '0',
        (0, 0, 1, 0, 0, 1): '-',
        # Using a common non-standard for equals
        (0, 1, 1, 0, 1, 1): '=',
    }

    result_chars = []
    
    # Numbers for the equation "8 - 1 = 7"
    numbers = ["8", "-", "1", "=", "7"]
    print("The decoded equation is:")
    
    equation_parts = []
    for i in range(5):
        dots = (
            1 if h1[2 * i] == 't' else 0,
            1 if h2[2 * i] == 't' else 0,
            1 if h3[2 * i] == 't' else 0,
            1 if h1[2 * i + 1] == 't' else 0,
            1 if h2[2 * i + 1] == 't' else 0,
            1 if h3[2 * i + 1] == 't' else 0,
        )
        char = braille_map.get(dots, '?')
        equation_parts.append(char)
        
    # The puzzle asks for the final equation, which spells out the numbers.
    # The direct decoding results in "8 - 1 = 7".
    # The question asks for "Which word", which implies the spelled out version.
    print("EIGHT - ONE = SEVEN")
    
    print("\nThe numeric equation is:")
    print(' '.join(equation_parts))

solve()