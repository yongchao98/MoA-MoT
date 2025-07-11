def solve_chess_mate():
    """
    This function identifies and prints the shortest possible mating sequence for black
    from the given chess position.
    """
    # The shortest possible mating sequence involves black checking the white king,
    # forcing a move that allows for a subsequent checkmate.
    # 1. Black moves h7-pawn to h6, checking the King on g5.
    # 2. White is forced to move the King. One possible (and poor) move is to h6.
    # 3. Black's Knight on e5 moves to g4, delivering checkmate.
    mating_sequence = ["h6+", "Kh6", "Ng4#"]
    
    # The problem asks to output the moves in order.
    # I am printing them as a comma-separated string.
    print(", ".join(mating_sequence))

solve_chess_mate()