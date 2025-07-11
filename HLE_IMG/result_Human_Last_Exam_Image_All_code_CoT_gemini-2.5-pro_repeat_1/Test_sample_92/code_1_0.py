def solve_art_puzzle():
    """
    Analyzes a list of paintings to determine which were made by a French painter before 1900.
    """
    paintings = [
        {'id': 'A', 'artist': 'Félix Vallotton', 'nationality': 'Swiss', 'year': 1892},
        {'id': 'B', 'artist': 'Édouard Vuillard', 'nationality': 'French', 'year': 1898},
        {'id': 'C', 'artist': 'Henri-Edmond Cross', 'nationality': 'French', 'year': 1895},
        {'id': 'D', 'artist': 'Odilon Redon', 'nationality': 'French', 'year': 1890},
        {'id': 'E', 'artist': 'Georges Rouault', 'nationality': 'French', 'year': 1906},
        {'id': 'F', 'artist': 'AI', 'nationality': 'N/A', 'year': 2022}, # Represents non-qualifying art
    ]
    
    result = ""
    for painting in paintings:
        is_french = painting['nationality'] == 'French'
        is_before_1900 = painting['year'] < 1900
        
        if is_french and is_before_1900:
            result += 'Y'
        else:
            result += 'N'
            
    print(result)

solve_art_puzzle()