def solve_art_puzzle():
    """
    Determines which artworks were created by a French painter before 1900.
    A: Félix Vallotton (Swiss, 1892) -> N
    B: Édouard Vuillard (French, c. 1908) -> N
    C: Henri-Edmond Cross (French, c. 1895) -> Y
    D: Odilon Redon (French, 1890) -> Y
    E: Henri de Toulouse-Lautrec (French, c. 1889) -> Y
    F: AI-generated (N/A, recent) -> N
    """
    
    # Each tuple contains: (Is French?, Is before 1900?)
    artwork_conditions = {
        'A': (False, True),  # Not French
        'B': (True, False),   # Not before 1900
        'C': (True, True),
        'D': (True, True),
        'E': (True, True),
        'F': (False, False)  # Not French and not before 1900
    }
    
    result = ""
    for letter in sorted(artwork_conditions.keys()):
        is_french, is_before_1900 = artwork_conditions[letter]
        if is_french and is_before_1900:
            result += "Y"
        else:
            result += "N"
            
    print(result)

solve_art_puzzle()