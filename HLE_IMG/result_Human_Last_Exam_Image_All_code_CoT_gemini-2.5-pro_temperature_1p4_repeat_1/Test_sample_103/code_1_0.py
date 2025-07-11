def calculate_material_advantage():
    """
    Analyzes the chess position after the winning combination for White.
    The sequence is 1. Nh6+ Kh8 2. Qxd8 Qxg2+ 3. Qxg2.
    """
    white_pieces = {
        "Queen": 9,
        "Rook": 5,
        "Knight": 3,
    }
    
    black_pieces = {
        "Rook": 5,
        "Bishop": 3,
    }

    white_material = sum(white_pieces.values())
    black_material = sum(black_pieces.values())
    
    advantage = white_material - black_material

    print("Analysis of the position after 1. Nh6+ Kh8 2. Qxd8 Qxg2+ 3. Qxg2:")
    print("-" * 30)
    print("White's remaining major pieces:")
    for piece, value in white_pieces.items():
        print(f"- {piece}: {value} points")
    print(f"Total for White: {white_material}")
    
    print("\nBlack's remaining major pieces:")
    for piece, value in black_pieces.items():
        print(f"- {piece}: {value} points")
    print(f"Total for Black: {black_material}")
    
    print("-" * 30)
    print("Final Material Advantage Calculation:")
    print(f"White's material ({white_material}) - Black's material ({black_material}) = {advantage} points")
    print("\nWhite has a decisive material advantage, leading to a winning position.")
    print("Therefore, the best move is Nh6+.")

calculate_material_advantage()