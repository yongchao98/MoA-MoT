def solve_art_puzzle():
    """
    Identifies which of the artworks A-F were created by a French painter before 1900.
    """
    paintings = {
        'A': {'title': 'Bathers on a Summer Evening', 'artist': 'Félix Vallotton', 'nationality': 'Swiss/French', 'year': 1892},
        'B': {'title': 'Homage to Cézanne (study)', 'artist': 'Maurice Denis', 'nationality': 'French', 'year': 1900},
        'C': {'title': 'The Olive Trees', 'artist': 'Henri-Edmond Cross', 'nationality': 'French', 'year': 1893},
        'D': {'title': 'Closed Eyes', 'artist': 'Odilon Redon', 'nationality': 'French', 'year': 1890},
        'E': {'title': 'La Buveuse (The Drinker)', 'artist': 'Émile Bernard', 'nationality': 'French', 'year': 1928},
        'F': {'title': 'Untitled Landscape', 'artist': 'Unknown/Modern Digital', 'nationality': 'N/A', 'year': 'c. 21st century'}
    }

    result_string = ""

    print("Analyzing the paintings based on the criteria: French painter AND created before 1900.\n")

    for key, data in sorted(paintings.items()):
        is_french = 'French' in data['nationality']
        
        # Check if year is a number and before 1900
        try:
            is_before_1900 = int(data['year']) < 1900
        except (ValueError, TypeError):
            is_before_1900 = False

        satisfies_conditions = is_french and is_before_1900
        
        if satisfies_conditions:
            result_char = 'Y'
            reason = "Yes (Artist is French, Year is before 1900)"
        else:
            result_char = 'N'
            if not is_french:
                reason = "No (Artist is not French)"
            elif not is_before_1900:
                reason = "No (Year is not before 1900)"

        print(f"Image {key}: '{data['title']}' by {data['artist']} ({data['nationality']}, {data['year']})")
        print(f"Result: {reason}\n")
        
        result_string += result_char

    print("Final answer format:")
    print(result_string)

solve_art_puzzle()