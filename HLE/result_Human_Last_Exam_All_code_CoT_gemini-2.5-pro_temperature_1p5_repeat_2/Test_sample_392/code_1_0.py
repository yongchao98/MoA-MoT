def solve_chess_puzzle():
    """
    This function prints the best sequence of moves for black to checkmate white
    from the given chess position.
    """
    black_move_1 = "1... Qg1+"
    white_move_2 = "2. Rxg1"
    black_move_2 = "2... Nf2#"

    print("The best sequence is a forced mate in 2, initiated by a queen sacrifice.")
    print("Here is the sequence of moves:")
    print(f"Black's first move: {black_move_1}")
    print(f"White's forced reply: {white_move_2}")
    print(f"Black's final move, delivering checkmate: {black_move_2}")
    print("\nThe full sequence in Algebraic Classic Notation is:")
    print(f"1... Qg1+ 2. Rxg1 Nf2#")

solve_chess_puzzle()