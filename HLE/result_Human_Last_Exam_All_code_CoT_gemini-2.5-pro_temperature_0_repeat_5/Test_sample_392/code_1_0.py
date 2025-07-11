def solve_chess_puzzle():
    """
    This function prints the best sequence of moves for black to checkmate white
    from the given chess position.
    """
    # The best sequence is C, which is a mate in 2 moves.
    # It's a famous pattern called Philidor's Legacy or a smothered mate.
    
    # The sequence is: 1... Qg1+  2. Rxg1  Nf2#
    
    # Define the moves in the sequence
    black_move_1 = "1... Qg1+"
    white_move_1 = "2. Rxg1"
    black_move_2 = "2... Nf2#"
    
    # Print the explanation and the sequence
    print("The best sequence is a forced mate in 2 for Black.")
    print("The moves are as follows:")
    print(f"Move 1 (Black): {black_move_1}")
    print(f"Move 2 (White): {white_move_1}")
    print(f"Move 2 (Black): {black_move_2}")
    
    # As requested, printing the final equation which represents the full sequence.
    # Each number (1, 2) from the notation is included.
    print("\nFinal Sequence:")
    print(f"{black_move_1} {white_move_1} {black_move_2}")

solve_chess_puzzle()