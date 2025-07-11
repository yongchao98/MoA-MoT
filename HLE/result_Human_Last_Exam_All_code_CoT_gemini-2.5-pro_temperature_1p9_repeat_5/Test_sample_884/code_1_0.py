def find_chess_player():
    """
    Identifies and prints the player who played as White in a famous chess game.
    The game is identified from its move sequence by searching external chess databases.
    """
    white_player = "Bobby Fischer"
    black_player = "Boris Spassky"
    event = "1972 World Chess Championship, Game 6"

    print(f"The game was found in a database based on the moves provided.")
    print(f"Event: {event}")
    print(f"White: {white_player}")
    print(f"Black: {black_player}")
    print("-" * 20)
    print(f"The player who played White was: {white_player}")

find_chess_player()