def solve_chess_puzzle():
    """
    This function prints the best sequence of moves for black to force a checkmate.
    The solution involves a classic smothered mate pattern.
    """
    moves = {
        "1...": "Nf2+",
        "2.": "Kg1",
        "2...": "Nh3+",
        "3.": "Kh1",
        "3...": "Qg1+",
        "4.": "Rxg1",
        "4...": "Nf2#"
    }
    
    print("The best sequence of moves is:")
    for turn, move in moves.items():
        print(f"{turn} {move}")

solve_chess_puzzle()
<<<E>>>