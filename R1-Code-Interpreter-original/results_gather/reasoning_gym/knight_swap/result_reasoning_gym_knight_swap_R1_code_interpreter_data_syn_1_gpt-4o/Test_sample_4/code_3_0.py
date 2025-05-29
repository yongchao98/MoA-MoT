from collections import deque

def is_valid_move(start, end):
    dx = abs(ord(end[0]) - ord(start[0]))
    dy = abs(int(end[1]) - int(start[1]))
    return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)

def swap_knights_bfs(initial_board):
    goal_state = {
        'w': ['A1', 'A2', 'A3'],
        'B': ['C1', 'C2']
    }
    
    queue = deque([(initial_board, [], 0)])
    visited = set()
    
    while queue:
        board, moves, turn = queue.popleft()
        
        # Check if the board is in the goal state
        if board['w'] == goal_state['w'] and board['B'] == goal_state['B']:
            return moves
        
        # Create a hashable state representation
        state = (tuple(board['w']), tuple(board['B']))
        if state in visited:
            continue
        visited.add(state)
        
        # Determine the current player
        current_player = 'w' if turn % 2 == 0 else 'B'
        
        # Try all possible moves for the current player
        for i, start in enumerate(board[current_player]):
            for end in ['A1', 'A2', 'A3', 'B1', 'B2', 'B3', 'C1', 'C2', 'C3']:
                if end not in board['w'] and end not in board['B'] and is_valid_move(start, end):
                    # Make the move
                    new_board = {k: v[:] for k, v in board.items()}
                    new_board[current_player][i] = end
                    # Add the new state to the queue
                    queue.append((new_board, moves + [f"{current_player},{start},{end}"], turn + 1))
    
    return None

# Initial board setup
initial_board = {
    'w': ['C1', 'C2'],
    'B': ['A1', 'A2', 'A3']
}

# Find the solution using BFS
solution = swap_knights_bfs(initial_board)

# Output the result
if solution:
    print(solution)
else:
    print("No")