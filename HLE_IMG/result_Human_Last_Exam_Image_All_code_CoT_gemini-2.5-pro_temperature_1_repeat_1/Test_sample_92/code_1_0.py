def solve_art_puzzle():
    """
    Identifies which of the six paintings were created by a French painter before 1900.
    """
    paintings = {
        'A': {'title': "Bathers on a Summer Evening", 'artist': "Félix Vallotton", 'nationality': "Swiss", 'year': 1892},
        'B': {'title': "Two Studies for 'A Modern Olympia'", 'artist': "Paul Cézanne", 'nationality': "French", 'year': 1874},
        'C': {'title': "The Olive Trees", 'artist': "Henri-Edmond Cross", 'nationality': "French", 'year': 1892},
        'D': {'title': "Closed Eyes", 'artist': "Odilon Redon", 'nationality': "French", 'year': 1890},
        'E': {'title': "Breton Woman with a Striped Blouse", 'artist': "Émile Bernard", 'nationality': "French", 'year': 1888},
        'F': {'title': "Modern Pointillism Landscape", 'artist': "Unknown/Modern", 'nationality': "Unknown", 'year': 2000} # Assuming post-2000
    }

    result_string = ""
    print("Analysis of each painting:")
    
    # Use sorted keys to ensure order A, B, C, D, E, F
    for label in sorted(paintings.keys()):
        painting = paintings[label]
        is_french = painting['nationality'] == "French"
        is_before_1900 = painting['year'] < 1900
        
        if is_french and is_before_1900:
            satisfies = True
            result_char = 'Y'
        else:
            satisfies = False
            result_char = 'N'
        
        result_string += result_char
        
        print(f"\nPainting {label}:")
        print(f"  - Artist: {painting['artist']} ({painting['nationality']})")
        print(f"  - Year: {painting['year']}")
        print(f"  - Is French? {'Yes' if is_french else 'No'}")
        print(f"  - Created before 1900? {'Yes' if is_before_1900 else 'No'}")
        print(f"  - Satisfies conditions? {'Yes' if satisfies else 'No'}")

    print("\nFinal Answer Format:")
    print(result_string)

solve_art_puzzle()