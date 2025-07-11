def solve_art_puzzle():
    """
    Analyzes six paintings to determine which were made by a French painter before 1900.
    """
    paintings = [
        {'id': 'A', 'title': 'Models', 'artist': 'Félix Vallotton', 'year': 1899, 'nationality': 'Swiss'},
        {'id': 'B', 'title': 'Chez le Père Martin', 'artist': 'Édouard Manet', 'year': 1874, 'nationality': 'French'},
        {'id': 'C', 'title': 'The Olive Trees', 'artist': 'Henri-Edmond Cross', 'year': 1892, 'nationality': 'French'},
        {'id': 'D', 'title': 'Closed Eyes', 'artist': 'Odilon Redon', 'year': 1890, 'nationality': 'French'},
        {'id': 'E', 'title': 'Old Woman', 'artist': 'Georges Rouault', 'year': 1905, 'nationality': 'French'},
        {'id': 'F', 'title': 'Modern Digital Art', 'artist': 'Unknown/AI', 'year': 2000, 'nationality': 'N/A'}
    ]

    result_string = ""
    print("Analyzing paintings...")
    print("-" * 30)

    for painting in paintings:
        is_french = painting['nationality'] == 'French'
        is_before_1900 = painting['year'] < 1900
        
        satisfies_conditions = is_french and is_before_1900
        
        result_char = 'Y' if satisfies_conditions else 'N'
        result_string += result_char
        
        print(f"Painting {painting['id']}:")
        print(f"  Artist: {painting['artist']} ({painting['nationality']})")
        print(f"  Year: {painting['year']}")
        print(f"  Condition Met (French Painter before 1900): {satisfies_conditions}")
        print("-" * 30)

    print(f"Final Answer String: {result_string}")

solve_art_puzzle()