def solve_art_puzzle():
    """
    Analyzes a list of artworks to determine which were made by French painters before 1900.
    """
    artworks = {
        'A': {'artist': 'Félix Vallotton', 'nationality': 'Swiss', 'year': 1899},
        'B': {'artist': 'Édouard Vuillard', 'nationality': 'French', 'year': 1898},
        'C': {'artist': 'Henri-Edmond Cross', 'nationality': 'French', 'year': 1892},
        'D': {'artist': 'Odilon Redon', 'nationality': 'French', 'year': 1890},
        'E': {'artist': 'Georges Rouault', 'nationality': 'French', 'year': 1907},
        'F': {'artist': 'Maximilien Luce', 'nationality': 'French', 'year': 1899}
    }

    result_string = ""
    # Sort the keys to ensure the order is A, B, C, D, E, F
    for label in sorted(artworks.keys()):
        artwork = artworks[label]
        is_french = artwork['nationality'] == 'French'
        is_before_1900 = artwork['year'] < 1900

        if is_french and is_before_1900:
            result_string += 'Y'
        else:
            result_string += 'N'
            
    print(result_string)

solve_art_puzzle()