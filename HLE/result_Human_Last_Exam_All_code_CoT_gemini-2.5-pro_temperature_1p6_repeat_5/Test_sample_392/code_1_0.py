def solve_chess_puzzle():
    """
    This function prints the best sequence of moves for Black to checkmate White
    from the given chess position.
    """
    
    # The best sequence is a classic "Smothered Mate" pattern.
    sequence = [
        "1... Qg1+",  # Black sacrifices the Queen. White must capture with the Rook.
        "2. Rxg1",    # White's forced move.
        "2... Nf2#"   # Black's Knight delivers checkmate.
    ]
    
    print("The best sequence of moves for Black to force checkmate is:")
    
    # We output each part of the final sequence, including the move numbers.
    print(f"Move 1 (Black): {sequence[0].split('... ')[1]}")
    print(f"Move 2 (White): {sequence[1].split('. ')[1]}")
    print(f"Move 2 (Black): {sequence[2].split('... ')[1]}")
    
    print("\nThe full sequence in Algebraic Notation is:")
    print(" ".join(sequence))

solve_chess_puzzle()