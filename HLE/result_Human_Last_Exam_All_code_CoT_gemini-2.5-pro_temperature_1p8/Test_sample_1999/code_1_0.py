def solve_spice_girls_riddle():
    """
    Determines which Spice Girl stands on the White Queen's starting square
    based on the riddle's conditions.
    """
    # Step 1: The order of the Spice Girls as referenced in the "Wannabe" rap bridge.
    # 1. "Em" -> Baby Spice
    # 2. "G" -> Ginger Spice
    # 3. "Easy V" -> Posh Spice
    # 4. "me" (rapped by Mel C) -> Sporty Spice
    # 5. The remaining member -> Scary Spice
    spice_girls_in_order = ["Baby", "Ginger", "Posh", "Sporty", "Scary"]

    # Step 2: The files (columns) on a chessboard.
    # The members are placed on rank 1, starting from the 'a' file.
    board_files = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']

    # The starting square for the White Queen.
    queen_start_square = "d1"

    print("Mapping the Spice Girls to their positions on the chessboard's first rank:")
    print("-" * 70)

    final_answer = ""
    # Step 3 & 4: Map each member to a square and find who is on the Queen's square.
    for i in range(len(spice_girls_in_order)):
        member = spice_girls_in_order[i]
        # Creates the square notation, e.g., 'a' + '1' -> 'a1'
        current_square = board_files[i] + '1'
        
        # This part fulfills the requirement to show the mapping for each member.
        print(f"The member at position {i+1} is {member} Spice, who stands on square {current_square}.")

        if current_square == queen_start_square:
            final_answer = member

    # Step 5: Announce the result.
    print("-" * 70)
    print(f"\nThe White Queen starts on square {queen_start_square}.")
    print(f"Therefore, the member standing on that square is {final_answer} Spice.")

solve_spice_girls_riddle()