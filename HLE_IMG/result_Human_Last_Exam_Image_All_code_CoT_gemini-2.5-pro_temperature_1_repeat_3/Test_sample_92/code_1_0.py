def solve_art_puzzle():
    """
    Analyzes a list of paintings to determine which were created by a French painter before 1900.
    """
    paintings = {
        'A': {'artist': 'Félix Vallotton', 'title': 'The Bathing', 'year': 1892, 'nationality_at_creation': 'Swiss'},
        'B': {'artist': 'Henri de Toulouse-Lautrec', 'title': 'In the Salon of the Rue des Moulins (study)', 'year': 1894, 'nationality_at_creation': 'French'},
        'C': {'artist': 'Henri-Edmond Cross', 'title': 'The Olive Trees', 'year': 1893, 'nationality_at_creation': 'French'},
        'D': {'artist': 'Odilon Redon', 'title': 'Closed Eyes', 'year': 1890, 'nationality_at_creation': 'French'},
        'E': {'artist': 'Georges Rouault', 'title': 'Prostitute', 'year': 1906, 'nationality_at_creation': 'French'},
        'F': {'artist': 'N/A (AI-generated)', 'title': 'Pointillist Landscape', 'year': 21, 'nationality_at_creation': 'N/A'}
    }

    result_string = ""
    for label in sorted(paintings.keys()):
        painting = paintings[label]
        is_french = painting['nationality_at_creation'] == 'French'
        is_before_1900 = painting['year'] < 1900
        
        if is_french and is_before_1900:
            result_string += 'Y'
        else:
            result_string += 'N'
            
    print(result_string)

solve_art_puzzle()