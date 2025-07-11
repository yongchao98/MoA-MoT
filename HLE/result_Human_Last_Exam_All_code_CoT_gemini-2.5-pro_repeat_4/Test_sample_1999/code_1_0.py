def solve_spice_girls_chess_riddle():
    """
    This function determines which Spice Girl would be on the White Queen's
    starting square based on the rap bridge from "Wannabe".
    """
    
    # 1. The order of the five Spice Girls as referenced in the song.
    # "Em" -> Emma (Baby)
    # "G" -> Geri (Ginger)
    # "V" -> Victoria (Posh)
    # "as for me" -> Mel B (Scary) is the one rapping this part.
    # The fifth member, to complete the group of five, is Mel C (Sporty).
    spice_girls_in_order = ["Baby", "Ginger", "Posh", "Scary", "Sporty"]
    
    # 2. The order of White's pieces on the first rank, starting from the queenside Rook.
    # The squares are a1, b1, c1, d1, e1...
    piece_positions = ["Queenside Rook", "Knight", "Bishop", "Queen", "King"]
    
    # 3. Map the Spice Girls to their positions on the chessboard.
    print("Placing the Spice Girls on the chessboard:")
    
    placement = {}
    for i in range(len(spice_girls_in_order)):
        member = spice_girls_in_order[i]
        position = piece_positions[i]
        placement[position] = member
        
        # This fulfills the "output each number in the final equation" requirement
        # by showing the mapping for each position.
        print(f"Position {i + 1} ({position}) = {member} Spice")
        
    # 4. Find which member is on the Queen's starting square.
    target_position = "Queen"
    member_on_queen_square = placement[target_position]
    
    print(f"\nThe White {target_position} starts at Position 4.")
    print(f"The Spice Girl at this position is {member_on_queen_square} Spice.")
    print("\nHer nickname prefix is the answer.")
    
    final_answer = member_on_queen_square
    print(f"\nFinal Answer: {final_answer}")

solve_spice_girls_chess_riddle()
<<<Scary>>>