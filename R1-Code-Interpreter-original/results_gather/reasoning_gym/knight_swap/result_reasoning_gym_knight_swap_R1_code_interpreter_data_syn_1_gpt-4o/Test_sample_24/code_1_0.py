def is_valid_move(start, end):
    dx = abs(ord(end[0]) - ord(start[0]))
    dy = abs(int(end[1]) - int(start[1]))
    return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)

def find_moves():
    # Initial positions
    white_knights = {'A3': 'w', 'D2': 'w'}
    black_knights = {'B3': 'B', 'C2': 'B'}
    empty_squares = {'A1', 'B1', 'C1', 'C3', 'B4', 'D4', 'E1', 'E4'}

    # Target positions
    target_white = {'B3', 'C2'}
    target_black = {'A3', 'D2'}

    # Moves list
    moves = []

    # Define possible moves for each knight
    knight_moves = {
        'A3': ['B1', 'C2'],
        'D2': ['B1', 'C3', 'E4'],
        'B3': ['A1', 'C1', 'D2'],
        'C2': ['A1', 'B4', 'D4', 'E1']
    }

    # Backtracking function to find the sequence of moves
    def backtrack(white_knights, black_knights, empty_squares, moves, depth):
        # Limit the depth to prevent infinite recursion
        if depth > 10:
            return None

        # Check if the swap is complete
        if set(white_knights.keys()) == target_white and set(black_knights.keys()) == target_black:
            return moves

        # Try moving each white knight
        for start in list(white_knights.keys()):
            for end in knight_moves[start]:
                if end in empty_squares:
                    # Make the move
                    moves.append(f"w,{start},{end}")
                    white_knights[end] = white_knights.pop(start)
                    empty_squares.remove(end)
                    empty_squares.add(start)

                    # Recurse
                    result = backtrack(white_knights, black_knights, empty_squares, moves, depth + 1)
                    if result:
                        return result

                    # Undo the move
                    moves.pop()
                    white_knights[start] = white_knights.pop(end)
                    empty_squares.add(end)
                    empty_squares.remove(start)

        # Try moving each black knight
        for start in list(black_knights.keys()):
            for end in knight_moves[start]:
                if end in empty_squares:
                    # Make the move
                    moves.append(f"B,{start},{end}")
                    black_knights[end] = black_knights.pop(start)
                    empty_squares.remove(end)
                    empty_squares.add(start)

                    # Recurse
                    result = backtrack(white_knights, black_knights, empty_squares, moves, depth + 1)
                    if result:
                        return result

                    # Undo the move
                    moves.pop()
                    black_knights[start] = black_knights.pop(end)
                    empty_squares.add(end)
                    empty_squares.remove(start)

        return None

    # Start the backtracking process
    result = backtrack(white_knights, black_knights, empty_squares, moves, 0)
    return result if result else "No"

# Execute the function and print the result
print(find_moves())