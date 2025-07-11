def solve_spice_girls_chess_puzzle():
    """
    Determines which Spice Girl stands on the White Queen's starting square
    based on the order from the "Wannabe" rap bridge.
    """

    # The order of the Spice Girls as referenced in the rap bridge.
    # 1. "Em" -> Emma Bunton (Baby Spice)
    # 2. "G" -> Geri Halliwell (Ginger Spice)
    # 3. "Easy V" -> Victoria Beckham (Posh Spice)
    # 4. "And as for me" -> Mel C, the singer of the line (Sporty Spice)
    # 5. Mel B (Scary Spice) is the rapper introducing the members and is conventionally placed last.
    spice_girls = ["Baby", "Ginger", "Posh", "Sporty", "Scary"]
    nicknames = [s + " Spice" for s in spice_girls]

    # The pieces on White's back rank, starting from the queenside (a-file).
    chess_pieces = ["Rook", "Knight", "Bishop", "Queen", "King"]
    files = ["a", "b", "c", "d", "e"]
    squares = [f + "1" for f in files]

    # The target piece we are looking for is the Queen.
    target_piece = "Queen"
    target_index = -1

    print("Mapping the Spice Girls to the chessboard squares:")
    print("-" * 50)
    for i in range(len(chess_pieces)):
        position = i + 1
        square = squares[i]
        piece = chess_pieces[i]
        spice_girl = nicknames[i]
        
        # This line corresponds to the "output each number" requirement.
        print(f"Position {position}: Square {square} ({piece}) -> {spice_girl}")
        
        if piece == target_piece:
            target_index = i

    print("-" * 50)
    
    if target_index != -1:
        result_nickname = spice_girls[target_index]
        print(f"The White Queen starts on square {squares[target_index]}.")
        print(f"The member standing on that square is {nicknames[target_index]}.")
        print("The answer is the word that precedes 'Spice' in her nickname.")
        print(f"\nFinal Answer: {result_nickname}")
    else:
        print("Could not determine the answer.")

solve_spice_girls_chess_puzzle()
