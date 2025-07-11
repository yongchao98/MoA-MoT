def solve_art_puzzle():
    """
    Analyzes a list of paintings to determine which were made by a French painter before 1900.
    """
    paintings = [
        {'id': 'A', 'artist': 'Félix Vallotton', 'nationality': 'Swiss', 'year': 1899, 'notes': 'Became French citizen in 1900.'},
        {'id': 'B', 'artist': 'Honoré Daumier', 'nationality': 'French', 'year': 1860},
        {'id': 'C', 'artist': 'Henri-Edmond Cross', 'nationality': 'French', 'year': 1895},
        {'id': 'D', 'artist': 'Odilon Redon', 'nationality': 'French', 'year': 1890},
        {'id': 'E', 'artist': 'Georges Rouault', 'nationality': 'French', 'year': 1906},
        {'id': 'F', 'artist': 'AI Generated', 'nationality': 'N/A', 'year': 2100, 'notes': 'Not a historical painting.'} # Year set to future to fail check
    ]

    result_string = ""
    for painting in paintings:
        # The conditions are: nationality must be 'French' and the year must be before 1900.
        is_french = painting['nationality'] == 'French'
        is_before_1900 = painting['year'] < 1900

        if is_french and is_before_1900:
            result_string += 'Y'
        else:
            result_string += 'N'

    print(result_string)

solve_art_puzzle()