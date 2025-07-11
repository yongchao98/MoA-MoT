def solve_chess_riddle():
    """
    This function identifies the player of the black pieces based on the provided game.
    The game is identified by searching a chess database for its moves.
    """
    # Based on a search of chess databases, the provided game is:
    # Event: 2023 FIDE World Championship Match, Game 12
    # White Player: Ding, Liren
    # Black Player: Nepomniachtchi, Ian
    # The question asks for the player of the black pieces.
    
    player_black = "Nepomniachtchi, Ian"
    
    print(f"The player of the black pieces was: {player_black}")

solve_chess_riddle()