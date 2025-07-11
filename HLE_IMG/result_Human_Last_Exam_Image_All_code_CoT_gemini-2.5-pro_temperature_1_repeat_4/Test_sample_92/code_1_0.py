def solve_art_puzzle():
    """
    Determines which paintings were created by a French painter before 1900.
    """
    paintings = [
        {'id': 'A', 'artist': 'Frédéric Bazille', 'nationality': 'French', 'year': 1869},
        {'id': 'B', 'artist': 'Pierre Bonnard', 'nationality': 'French', 'year': 1904},
        {'id': 'C', 'artist': 'Henri-Edmond Cross', 'nationality': 'French', 'year': 1897},
        {'id': 'D', 'artist': 'Odilon Redon', 'nationality': 'French', 'year': 1890},
        {'id': 'E', 'artist': 'Georges Rouault', 'nationality': 'French', 'year': 1910},
        {'id': 'F', 'artist': 'Georges Seurat', 'nationality': 'French', 'year': 1888}
    ]

    result_string = ""
    print("Checking which paintings were created by a French painter before 1900:\n")

    for painting in paintings:
        is_french = painting['nationality'] == 'French'
        is_before_1900 = painting['year'] < 1900
        
        passes = is_french and is_before_1900
        result_char = 'Y' if passes else 'N'
        result_string += result_char
        
        print(f"Painting {painting['id']}:")
        print(f"  Artist: {painting['artist']} ({painting['nationality']}), Year: {painting['year']}")
        print(f"  Is French? {is_french}")
        print(f"  Is before 1900? {is_before_1900}")
        print(f"  Result: {result_char}\n")

    print("Final Answer String:")
    print(result_string)

solve_art_puzzle()