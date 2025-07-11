def find_minimal_mate_configuration():
    """
    Calculates and prints the minimal piece configuration for a specific
    Diagonal Corridor Mate, adhering to all problem constraints.

    The logic, as explained above, concludes that a 3-piece solution is
    the minimum that can satisfy the "Black has material advantage" rule.
    This solution consists of:
    - White: A Bishop (value 3) as the attacker.
    - Black: A Pawn (value 1) and a Knight (value 3) as blockers.
    This configuration results in a material score of White 3 vs. Black 4.
    """
    
    # Define piece values for sorting
    piece_values = {
        "Pawn": 1,
        "Knight": 3,
        "Bishop": 3,
        "Rook": 5,
        "Queen": 9
    }

    # The pieces determined by the logic
    white_pieces = ["Bishop"]
    black_pieces = ["Pawn", "Knight"]

    # Sort each list by piece value
    white_pieces.sort(key=lambda p: piece_values[p])
    black_pieces.sort(key=lambda p: piece_values[p])

    # Format the final list of piece names
    final_piece_list = []
    for piece in white_pieces:
        final_piece_list.append(f"White {piece}")
    for piece in black_pieces:
        final_piece_list.append(f"Black {piece}")

    # Join the list into a comma-separated string and print
    result = ", ".join(final_piece_list)
    print(result)

find_minimal_mate_configuration()
<<<White Bishop, Black Pawn, Black Knight>>>