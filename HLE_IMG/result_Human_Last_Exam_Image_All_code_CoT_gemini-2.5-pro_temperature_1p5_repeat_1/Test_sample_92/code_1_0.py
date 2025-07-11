def solve_art_puzzle():
    """
    Analyzes a list of artworks to determine which were created by French painters before 1900.
    """
    artworks = [
        {'id': 'A', 'artist': 'Félix Vallotton', 'nationality': 'Swiss', 'year': 1897},
        {'id': 'B', 'artist': 'Edgar Degas', 'nationality': 'French', 'year': 1879},
        {'id': 'C', 'artist': 'Henri-Edmond Cross', 'nationality': 'French', 'year': 1892},
        {'id': 'D', 'artist': 'Odilon Redon', 'nationality': 'French', 'year': 1890},
        {'id': 'E', 'artist': 'Émile Bernard', 'nationality': 'French', 'year': 1888},
        {'id': 'F', 'artist': 'AI-generated', 'nationality': 'N/A', 'year': 2020} # Approximate year
    ]

    result = ""
    for art in artworks:
        # Check if the painter is French AND the year is before 1900
        if art['nationality'] == 'French' and art['year'] < 1900:
            result += 'Y'
        else:
            result += 'N'
    
    print(result)

solve_art_puzzle()