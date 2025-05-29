def is_valid_move(start, end):
    # Calculate the difference in positions
    dx = abs(ord(end[0]) - ord(start[0]))
    dy = abs(int(end[1]) - int(start[1]))
    # Check if the move is a valid knight move
    return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)

def generate_moves():
    # Initial positions
    white_knights = {'C3', 'D1'}
    black_knights = {'A2', 'B1'}
    empty_squares = {'.', 'B2'}

    # Possible moves
    moves = []

    # Try to find a sequence of moves
    # White moves first
    for w_start in white_knights:
        for w_end in empty_squares:
            if is_valid_move(w_start, w_end):
                # Move white knight
                moves.append(f"w,{w_start},{w_end}")
                white_knights.remove(w_start)
                white_knights.add(w_end)
                empty_squares.remove(w_end)
                empty_squares.add(w_start)

                # Black moves
                for b_start in black_knights:
                    for b_end in empty_squares:
                        if is_valid_move(b_start, b_end):
                            # Move black knight
                            moves.append(f"B,{b_start},{b_end}")
                            black_knights.remove(b_start)
                            black_knights.add(b_end)
                            empty_squares.remove(b_end)
                            empty_squares.add(b_start)

                            # Check if the swap is complete
                            if white_knights == {'A2', 'B1'} and black_knights == {'C3', 'D1'}:
                                return moves

                            # Revert black move
                            moves.pop()
                            black_knights.remove(b_end)
                            black_knights.add(b_start)
                            empty_squares.remove(b_start)
                            empty_squares.add(b_end)

                # Revert white move
                moves.pop()
                white_knights.remove(w_end)
                white_knights.add(w_start)
                empty_squares.remove(w_start)
                empty_squares.add(w_end)

    return "No"

# Generate and print the moves
print(generate_moves())