def solve_art_puzzle():
    """
    Analyzes a list of paintings to determine which were made by a French painter before 1900.
    """
    paintings = {
        'A': {'artist': 'Georges Seurat', 'nationality': 'French', 'year': 1888},
        'B': {'artist': 'Pierre Bonnard', 'nationality': 'French', 'year': 1905},
        'C': {'artist': 'Henri-Edmond Cross', 'nationality': 'French', 'year': 1892},
        'D': {'artist': 'Odilon Redon', 'nationality': 'French', 'year': 1890},
        'E': {'artist': 'Henri de Toulouse-Lautrec', 'nationality': 'French', 'year': 1889},
        'F': {'artist': 'AI Generator', 'nationality': 'N/A', 'year': 2022},
    }

    result = ""
    for painting_id in sorted(paintings.keys()):
        details = paintings[painting_id]
        # Check if the painter is French and the year is before 1900
        if details['nationality'] == 'French' and details['year'] < 1900:
            result += 'Y'
        else:
            result += 'N'
            
    print(result)

solve_art_puzzle()