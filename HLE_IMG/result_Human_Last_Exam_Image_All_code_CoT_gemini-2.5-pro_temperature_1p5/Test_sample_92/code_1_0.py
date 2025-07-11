def solve_painting_query():
    """
    Analyzes six paintings to determine which were created by a French painter before 1900.
    """
    paintings_data = {
        'A': {'title': 'The Models', 'artist': 'Félix Vallotton', 'nationality': 'Swiss', 'year': 1897},
        'B': {'title': 'Study for Homage to Cézanne', 'artist': 'Maurice Denis', 'nationality': 'French', 'year': 1900},
        'C': {'title': 'The Olive Trees', 'artist': 'Henri-Edmond Cross', 'nationality': 'French', 'year': 1892},
        'D': {'title': 'Closed Eyes', 'artist': 'Odilon Redon', 'nationality': 'French', 'year': 1890},
        'E': {'title': 'The Hangover', 'artist': 'Henri de Toulouse-Lautrec', 'nationality': 'French', 'year': 1889},
        'F': {'title': 'Modern Digital Art', 'artist': 'Unknown', 'nationality': 'N/A', 'year': '21st Century'}
    }

    result = ""
    for label in sorted(paintings_data.keys()):
        painting = paintings_data[label]
        is_french = (painting['nationality'] == 'French')
        
        is_before_1900 = False
        # Handle cases where year is a single integer or a string like '21st Century'
        if isinstance(painting['year'], int):
            is_before_1900 = painting['year'] < 1900

        if is_french and is_before_1900:
            result += 'Y'
        else:
            result += 'N'
    
    print(result)

solve_painting_query()