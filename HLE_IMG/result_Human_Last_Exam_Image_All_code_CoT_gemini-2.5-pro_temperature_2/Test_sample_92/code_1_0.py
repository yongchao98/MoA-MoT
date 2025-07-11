def solve_art_puzzle():
    """
    Analyzes six paintings to determine if they were made by a French painter before 1900.

    The function stores the details of each painting and then evaluates them against the given criteria.
    It prints the final result as a six-character string of 'Y's and 'N's.
    """
    paintings = {
        'A': {'artist': 'Félix Vallotton', 'nationality': 'Swiss/French', 'year': 1899, 'info': 'Active in France, part of French Les Nabis group'},
        'B': {'artist': 'Jean-Louis Forain', 'nationality': 'French', 'year': 1885, 'info': 'Approximate date for sketch'},
        'C': {'artist': 'Henri-Edmond Cross', 'nationality': 'French', 'year': 1893, 'info': ''},
        'D': {'artist': 'Odilon Redon', 'nationality': 'French', 'year': 1890, 'info': 'Date of the original famous version'},
        'E': {'artist': 'Émile Bernard', 'nationality': 'French', 'year': 1887, 'info': ''},
        'F': {'artist': 'AI Generated', 'nationality': 'N/A', 'year': 2023, 'info': 'Not a historical painting'},
    }

    result = []
    for key in sorted(paintings.keys()):
        painting = paintings[key]
        
        # Condition 1: Was the painter French (or primarily active in France)?
        is_french = 'French' in painting['nationality']
        
        # Condition 2: Was the painting created before 1900?
        is_before_1900 = painting['year'] < 1900
        
        if is_french and is_before_1900:
            result.append('Y')
        else:
            result.append('N')
            
    print("".join(result))

solve_art_puzzle()