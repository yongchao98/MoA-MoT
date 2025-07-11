def solve_art_puzzle():
    """
    Analyzes a list of artworks to determine which were created by a French painter before 1900.
    """
    artworks = [
        {'id': 'A', 'artist': 'Félix Vallotton', 'nationality': 'Swiss', 'year': 1892},
        {'id': 'B', 'artist': 'Édouard Manet', 'nationality': 'French', 'year': 1881},
        {'id': 'C', 'artist': 'Henri-Edmond Cross', 'nationality': 'French', 'year': 1895},
        {'id': 'D', 'artist': 'Odilon Redon', 'nationality': 'French', 'year': 1890},
        {'id': 'E', 'artist': 'Henri de Toulouse-Lautrec', 'nationality': 'French', 'year': 1889},
        {'id': 'F', 'artist': 'AI / Digital Artist', 'nationality': 'N/A', 'year': 2020} # Year is a placeholder for modern art
    ]

    result_string = ""
    for art in artworks:
        # Condition 1: Was the painter French?
        # Note: For artist A, Félix Vallotton became a French citizen in 1900, but was Swiss when he painted the piece in 1892.
        is_french = art['nationality'] == 'French'
        
        # Condition 2: Was the painting created before the year 1900?
        is_before_1900 = art['year'] < 1900

        if is_french and is_before_1900:
            result_string += 'Y'
        else:
            result_string += 'N'
            
    print(result_string)

solve_art_puzzle()