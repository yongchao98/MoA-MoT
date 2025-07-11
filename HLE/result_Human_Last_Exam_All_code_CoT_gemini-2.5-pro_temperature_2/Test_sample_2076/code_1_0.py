def solve_chess_puzzle():
    """
    This function determines and prints the minimal chess position for a Diagonal Corridor Mate
    based on the specified constraints.
    """

    # The pieces determined by the logic outlined in the explanation.
    # Format: (Name, Value)
    white_pieces = [("Bishop", 3)]
    black_pieces = [("Pawn", 1), ("Knight", 3)]

    # Sort each list by piece value in ascending order.
    # While these lists are already sorted, this step ensures correctness.
    white_pieces.sort(key=lambda p: p[1])
    black_pieces.sort(key=lambda p: p[1])

    # Construct the list of piece descriptions for the final output.
    output_parts = []
    for piece_name, _ in white_pieces:
        output_parts.append(f"White {piece_name}")

    for piece_name, _ in black_pieces:
        output_parts.append(f"Black {piece_name}")

    # Join the parts into a single comma-separated string.
    final_answer_string = ", ".join(output_parts)

    print(final_answer_string)

solve_chess_puzzle()