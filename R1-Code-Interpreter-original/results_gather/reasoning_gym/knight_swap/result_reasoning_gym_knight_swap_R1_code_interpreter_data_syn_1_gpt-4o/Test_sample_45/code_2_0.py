def is_valid_move(start, end):
    # Calculate the difference in positions
    dx = abs(ord(end[0]) - ord(start[0]))
    dy = abs(int(end[1]) - int(start[1]))
    # Check if the move is a valid knight move
    return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)

def find_moves():
    # Initial positions
    white_knights = ['B1', 'D1']
    black_knights = ['C3', 'B2']
    empty_squares = ['A3', 'C2', 'D3']

    # Target positions
    target_white = ['C3', 'B2']
    target_black = ['B1', 'D1']

    # Possible moves
    moves = []

    # Set to track visited states
    visited = set()

    # Try to find a sequence of moves
    def backtrack(w_pos, b_pos, empty, move_list, turn):
        # Sort positions to create a unique state representation
        state = (tuple(sorted(w_pos)), tuple(sorted(b_pos)), tuple(sorted(empty)))
        
        # Base case: check if the current positions match the target positions
        if sorted(w_pos) == sorted(target_white) and sorted(b_pos) == sorted(target_black):
            return move_list

        # Check if the state has been visited
        if state in visited:
            return None
        visited.add(state)

        if turn == 'w':
            for i, pos in enumerate(w_pos):
                for e in empty:
                    if is_valid_move(pos, e):
                        new_w_pos = w_pos[:i] + [e] + w_pos[i+1:]
                        new_empty = empty[:]
                        new_empty.remove(e)
                        new_empty.append(pos)
                        result = backtrack(new_w_pos, b_pos, new_empty, move_list + [f"w,{pos},{e}"], 'b')
                        if result:
                            return result
        else:
            for i, pos in enumerate(b_pos):
                for e in empty:
                    if is_valid_move(pos, e):
                        new_b_pos = b_pos[:i] + [e] + b_pos[i+1:]
                        new_empty = empty[:]
                        new_empty.remove(e)
                        new_empty.append(pos)
                        result = backtrack(w_pos, new_b_pos, new_empty, move_list + [f"B,{pos},{e}"], 'w')
                        if result:
                            return result
        return None

    # Start the backtracking with white's turn
    result = backtrack(white_knights, black_knights, empty_squares, moves, 'w')
    return result

# Find the sequence of moves
solution = find_moves()
print(solution if solution else "No")