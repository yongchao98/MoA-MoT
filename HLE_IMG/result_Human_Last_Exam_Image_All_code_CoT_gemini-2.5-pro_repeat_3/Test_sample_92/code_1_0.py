def solve_art_puzzle():
    """
    Analyzes a list of artworks to determine which were made by French painters before 1900.
    """
    artworks = [
        {'id': 'A', 'artist': 'Félix Vallotton', 'nationality': 'Swiss-French', 'year': 1897},
        {'id': 'B', 'artist': 'Édouard Vuillard', 'nationality': 'French', 'year': 1895},
        {'id': 'C', 'artist': 'Henri-Edmond Cross', 'nationality': 'French', 'year': 1892},
        {'id': 'D', 'artist': 'Odilon Redon', 'nationality': 'French', 'year': 1890},
        {'id': 'E', 'artist': 'Georges Rouault', 'nationality': 'French', 'year': 1906},
        {'id': 'F', 'artist': 'Modern/AI', 'nationality': 'N/A', 'year': 2023}
    ]

    result = []
    for art in artworks:
        # A painter is considered French if their nationality contains "French".
        is_french = 'French' in art['nationality']
        # The work must be created before the year 1900.
        is_before_1900 = art['year'] < 1900

        if is_french and is_before_1900:
            result.append('Y')
        else:
            result.append('N')

    print("".join(result))

solve_art_puzzle()