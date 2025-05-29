def is_valid_move(start, end):
    # Calculate the move difference
    dx = abs(ord(end[0]) - ord(start[0]))
    dy = abs(int(end[1]) - int(start[1]))
    # Check if the move is a valid knight move
    return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)

def swap_knights(board, moves, turn, visited):
    # Check if the board is in the goal state
    if board['w'] == ['A1', 'A2', 'A3'] and board['B'] == ['C1', 'C2']:
        return moves

    # Create a hashable state representation
    state = (tuple(board['w']), tuple(board['B']))
    if state in visited:
        return None
    visited.add(state)

    # Determine the current player
    current_player = 'w' if turn % 2 == 0 else 'B'

    # Try all possible moves for the current player
    for i, start in enumerate(board[current_player]):
        for end in ['A1', 'A2', 'A3', 'B1', 'B2', 'B3', 'C1', 'C2', 'C3']:
            if end not in board['w'] and end not in board['B'] and is_valid_move(start, end):
                # Make the move
                board[current_player][i] = end
                # Recurse with the new board state
                result = swap_knights(board, moves + [f"{current_player},{start},{end}"], turn + 1, visited)
                if result:
                    return result
                # Undo the move
                board[current_player][i] = start

    visited.remove(state)
    return None

# Initial board setup
initial_board = {
    'w': ['C1', 'C2'],
    'B': ['A1', 'A2', 'A3']
}

# Find the solution
solution = swap_knights(initial_board, [], 0, set())

# Output the result
if solution:
    print(solution)
else:
    print("No")